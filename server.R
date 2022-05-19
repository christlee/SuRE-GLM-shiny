library(rtracklayer)
library(data.table)
library(shiny)
library(AnnotationDbi)
library(ggplot2)
library(cowplot)
library(grid)
library(stringr)
library(here)
library(DT)
library(shinyjs)

peak_file  = c('K562' = 'data/SuRE_elasticNet_peaks_K562.bed',
               'HEPG2' = 'data/SuRE_elasticNet_peaks_K562.bed')

tss_file = 'data/tss_matrix.txt.gz'
lookup_matrix = 'data/lookup_matrix.txt.gz'


offset_coefs = readRDS("data/optimalFitOffsetCoefficients_K562_spatial.rds")


ucsc = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?"

ucsc_opt = paste("db=hg19",
                 "hubUrl=https://raw.githubusercontent.com/christlee/SuRE-GLM-shiny/main/trackhub.txt",
                 sep='&')

hg38ToHg19 = "data/hg38ToHg19.over.chain"
hg19ToHg38 = "data/hg19ToHg38.over.chain"


tss_gr = makeGRangesFromDataFrame(fread(cmd=paste("zcat", tss_file)),
                                  keep.extra.columns=T)
names(tss_gr) = tss_gr$tx_name

lookup_dt = fread(cmd=paste("zcat", lookup_matrix))
lookup_dt[, gene_name:=toupper(gene_name)]
lookup_dt[, ensembl_id := gsub('[.].*', '', gene_id)]

ensembl_dt = lookup_dt[,list(ensembl_id, transcript_id, gene_name)]
setkey(ensembl_dt, 'ensembl_id')


gencode_dt = lookup_dt[,list(gene_id = gsub('_.*', '', gene_id),
                             transcript_id)]
setkey(gencode_dt, 'gene_id')

symbol_dt = lookup_dt[, c('gene_name', 'transcript_id', 'ensembl_id')]
setkey(symbol_dt, 'gene_name')

transcript_dt = data.table(ensembl_tid = gsub('[.]*', '', names(tss_gr)),
                           transcript_id = names(tss_gr), key='ensembl_tid',
                           stringsAsFactors=F)


lookup_list = list('symbol'=symbol_dt, 'ensembl_gene'=ensembl_dt,
                   'ensembl_transcript'=transcript_dt, 'gencode'=gencode_dt)


get_center <- function(text_input, lookup_list, tss_gr){
  input_upper = toupper(text_input)
  match_list = list(ensembl_transcript='^ENST[0-9]+$',
                    ensembl_gene='^ENSG[0-9]+$',
                    gencode='^ENSG[0-9]+[.][0-9]+$',
                    symbol='.*')
  if (grepl('chr[0-9XY]+:[0-9]+:[-+]$', text_input)){
    split_vec = strsplit(text_input, ':')[[1]]
    center = GRanges(seqnames=split_vec[1],
                     IRanges(start=as.numeric(split_vec[2]),
                             width=1),
                     strand = split_vec[3])
  } else if (grepl('chr[0-9XY]+:[0-9]+$', text_input)){
      split_vec = strsplit(text_input, ':')[[1]]
      center = GRanges(seqnames=split_vec[1],
                       IRanges(start=as.numeric(split_vec[2]),
                               width=1),
                       strand='+')
  } else if (grepl('ENST[0-9]+[.][0-9]+', input_upper)){
    center = tss_gr[input_upper]
  } else {
    i = 1
    is_match = F
    while (!is_match==T) {
      is_match = grepl(match_list[i], input_upper)
      match_name = names(match_list)[i]
      i = i + 1
    }
    transcript_id = lookup_list[[match_name]][input_upper, transcript_id]

    if (!is.na(transcript_id)){
        center = tss_gr[transcript_id]
    } else {
        center = NA
    }

  }
  return(center)
}

#
triangle_dt <- function(center, upstream = 1000, downstream=1000,
                        cutoff = 2000, binsize=10, offset=0,
                        path = "data/",
                        filePart = "SURE_elasticNet_allele_K562_",
                        rev_comp = F) {
    start_time <- Sys.time()
    center_pos <- start(center)
    xregion <- promoters(center, upstream=upstream,
                         downstream=downstream)
    if (upstream + downstream < cutoff){
      region <- resize(xregion, width=cutoff, fix="center")
    } else {
      region <- xregion
    }

    half_bin = round(binsize/2)
    start_i = (end(region) - start(region) - half_bin) %% binsize
    if (start_i < half_bin){
        start(region) = start(region) + start_i - half_bin
    }
    if (rev_comp==F){
        strand <- ifelse(as.character(strand(center)) == "+", "plus", "minus")
    } else {
        strand <- ifelse(as.character(strand(center)) == "+", "minus", "plus")
    }
    file <- paste0(path, filePart, strand, ".bw")
    suppressWarnings(glm <- import(file, selection = BigWigSelection(region),
                                   as = "NumericList")[[1]] + offset)

    i <- seq(((end(region) - start(region) - half_bin) %% binsize) + 1,
             end(region) - start(region), binsize)

    xlim <- c(start(xregion), end(xregion))

    if(strand == "minus"){
        i <- rev(i)
        xlim <- rev(xlim)
    }
    i_list = lapply(i, function(n){(n-half_bin):(n+half_bin-1)})
    ## to make use of fast matrix operations, we first are going to
    ## create a slightly tilted triangle in a matrix.
    ## if we have a triangle like this:

    ##         10
    ##       8   9
    ##     5   6   7
    ##   1   2   3   4

    ## in the matrix, it will look like this
    ## 1  5  8  10
    ##    2  6  9
    ##       3  7
    ##          4

    max_bin_i = round(cutoff / binsize)
    mat <- matrix(0, length(i), length(i))
    glm_vec = vector('numeric', length=length(i))

    for(n in 1:length(i)){
        glm_vec[n] = sum(glm[i_list[[n]]])
        mat[1:n,n:length(i)] <- mat[1:n,n:length(i)] + sum(glm[i_list[[n]]])
        mat[n,1:(n-1)] <- NA
        mat[n,(1:length(i)) > n + max_bin_i] <- NA
    }

    start_mat <- matrix(start(region) + i - 1, nrow = length(i), ncol = length(i))

    end_mat <- matrix(start(region) + i - 1, nrow = length(i), ncol = length(i), byrow = T)

    end_time <- Sys.time()
    print(end_time - start_time)
    x_mat <- (start_mat + end_mat)[!is.na(mat)]/2
    y_mat <- abs(end_mat - start_mat)[!is.na(mat)] + binsize
    #
    # dt = data.table(x=x_mat, y=y_mat, score=mat[!is.na(mat)],
    #                 x_rel = (x_mat - center_pos) * ifelse(strand=='plus', -1, 1))
    #


    x_vec = (start(region):end(region))[i]
    # flat_dt = data.table(score=glm_vec,
    #                      x = x_vec,
    #                      x_rel = (x_vec - center_pos) * ifelse(strand=='plus', -1, 1))
    #
    ##TODO: make sure the spacing up and downstream don't get considered in the
    ##      calculation
    if ((strand=='plus' & rev_comp==F) | (strand=='minus' & rev_comp==T)){
        dt = data.table(x=x_mat, y=y_mat, score=mat[!is.na(mat)],
                        x_rel = (x_mat - center_pos))
        mat_dt = dt[x_rel <= downstream & x_rel >= -upstream, ]

        flat_dt = data.table(score=glm,
                             x = start(region):end(region))
        flat_dt[,x_rel := x - center_pos]
    } else {
        dt = data.table(x=x_mat, y=y_mat, score=mat[!is.na(mat)],
                        x_rel = (center_pos - x_mat))
        mat_dt = dt[x_rel <= downstream & x_rel >= -upstream, ]


        flat_dt = data.table(score=glm,
                             x = start(region):end(region))
        flat_dt[, x_rel := center_pos - x]
    }
    setkey(flat_dt, 'x')
    return(list(mat_dt=mat_dt,
                flat_dt=flat_dt,
                region=region))
}

lib_vec = c('K562'='SURE_elasticNet_allele_K562_',
            'HEPG2'='SURE_elasticNet_allele_HEPG2_')





input = list(ROI='NUP214', binsize=10, window=c(-1000,1000), cutoff = c(200,1500))

# colorlut <- colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
#                                "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61",
#                                "#F46D43", "#D73027", "#B91326",
#                                "#A50026"))(101) # height color lookup table
colorlut <- c(colorRampPalette(c("#ffffff", "#000000"))(91), rep("#000000",10)) # height color lookup table


plot_mat <-function(mat, input, cutoff, max_color, binsize=10){
  min_rel = mat[, min(x_rel)]
  max_rel = mat[, max(x_rel)]
  max_y = mat[,max(y)]
  corners = data.frame(x_rel=c(min_rel, min_rel + cutoff/2,
                               max_rel-max_y/2, max_rel),
                       y=c(0,max_y, max_y, 0))
  trapezoid = cbind(corners[c(1,2,2,3,3,4,4,1), ], line_n=c(1,1,2,2,3,3,4,4))
  p = ggplot(mat, aes(x=x_rel, y=y, fill=exp(score))) +
    geom_tile(width=binsize, height=binsize) +
    geom_line(data=trapezoid, aes(x=x_rel,y=y,group=line_n), inherit.aes=F) +
    ggplot2::theme_bw() +
    scale_fill_gradientn(colours=colorlut, limits=c(0,exp(max_color))) +
    ylab('reporter length') +
    xlab('##chromosome##') +
    coord_fixed(ylim=c(0,cutoff)) +
    ggplot2::theme(axis.title.x=element_blank(),
                   axis.text.y=element_text(size=14),
                   axis.title.y=element_text(size=14),
                   axis.text.x=element_text(size=14),
                   plot.margin=margin(t=0,r=20,b=0,l=0),
                   plot.title = element_text(hjust=0.5, size=14))
  return(p)
}

grTodt <- function(gr, center){
    dt = as.data.table(gr)
    dt_melt = melt(dt, measure.vars=c('start', 'end'), value.name='x')
    dt_melt[,xrel := x - center]
    dt_melt[,tx_name:=unlist(tx_name)]
    return(data.table(dt_melt))
}


printRegion <- function(region){
    paste0(region[1], ':', region[2], '-', region[3])
}

liftOver <- function(region, toHg19=TRUE){
    chain = ifelse(toHg19, hg38ToHg19, hg19ToHg38)
    region_str = paste(region, collapse='\t')

    cmd = paste0('printf "', region_str, '" | ', "liftOver",
                 " /dev/stdin ", here(chain), " /dev/stdout ", tempfile())
    print(cmd)
    lift_str = system(cmd, intern=T, ignore.stderr=T)
    lift_vec = strsplit(lift_str, '\t')[[1]]
    return(lift_vec)
}


get_peaks <- function(region, peak_file){
    region_str = paste(seqnames(region), start(region), end(region), sep='\t')
    cmd = paste0("printf '", region_str, "' | ",
                 "bedtools intersect -wa ",
                 "-a ", peak_file, " -b -")

    return(fread(cmd=cmd, col.names=c('seqnames', 'start', 'end', 'name',
                                      'score', 'strand')))
}





shinyServer(function(input, output, session) {
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
    updateSelection <- function(){
        updateTextInput(session, "hg19_selection",
            value=printRegion(vals$region))
        updateTextInput(session, "hg38_selection",
            value=printRegion(liftOver(vals$region, toHg19=F)))

    }

    vals <- reactiveValues(
        x = NA,
        y = NA,
        start = NA,
        end = NA,
        region = NA,
        center = NA,
        binsize = 10,
        show_minus = FALSE,
        uptodate = TRUE,
        fragment_warning = FALSE
    )


    observeEvent(input$window, {
        vals$uptodate = FALSE
    })

    observeEvent(input$lib, {
        vals$uptodate = FALSE
    })

    observeEvent(input$show_minus, {
        vals$uptodate = FALSE
    })


    dataInput <- eventReactive(input$go, {
        vals$uptodate = TRUE
        center = get_center(input$ROI, lookup_list, tss_gr)
        validate(
            need(!is.na(center), "please enter a valid gene name/identifier")
        )
        vals$center = center
        if (input$lib%in%c('K562', 'HEPG2')){
          offset = offset_coefs$width[[1]]
          cutoff = input$cutoff42
        } else {
          offset = 0
          cutoff = input$cutoff23
        }

        dt = triangle_dt(center, upstream=input$window[1]*-1,
                         downstream=input$window[2],
                         cutoff = cutoff, binsize=10, offset=offset,
                         filePart = lib_vec[input$lib])
        max_sense = dt$mat_dt[, max(score)]

        if (input$show_minus==TRUE){
            vals$show_minus = TRUE
            dt_rev = triangle_dt(center, upstream=input$window[1]*-1,
                                 downstream=input$window[2],
                                 cutoff = cutoff, binsize=vals$binsize, offset=offset,
                                 filePart = lib_vec[input$lib], rev_comp = T)
            names(dt_rev) = paste0(names(dt_rev), '_rev')
            max_antisense = dt_rev$mat_dt_rev[, max(score)]
        } else {
            vals$show_minus = FALSE
            dt_rev = NULL
            max_antisense = -Inf
        }
        strand = strand(center)
        sec = start(center) * ifelse(strand=='+', 1, -1)


        # peak_dt = suppressWarnings(get_peaks(dt$region,
        #                                      peak_file[input$lib]))
        # if (nrow(peak_dt)==0){
        #     peak_melt = peak_dt
        # } else {
        #     peak_dt[,group:=1:nrow(peak_dt)]
        #     peak_melt = melt(peak_dt, measure.vars=c('start', 'end'),
        #                      value.name='x')
        #     if (as.character(strand) == '+'){
        #         peak_melt[,x_rel := x - start(center)]
        #         peak_dt[, strand_order:= ifelse(strand=='+',2,1)]
        #         max_peak = peak_dt[order(strand_order,score, decreasing=T), ][1,]
        #         updateVals(start=max_peak[['start']] - start(center),
        #                    end=max_peak[['end']] - start(center))
        #     } else {
        #         peak_melt[,x_rel := start(center) - x]
        #         peak_dt[, strand_order:= ifelse(strand=='+',1,2)]
        #         max_peak = peak_dt[order(strand_order,score, decreasing=T), ][1,]
        #
        #         updateVals(start=start(center) - max_peak[['start']],
        #                    end=start(center) - max_peak[['end']])
        #     }
        # }
        max_color = max(max_sense, max_antisense)
        p = plot_mat(dt$mat_dt, input, cutoff, max_color, vals$binsize) +
          scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0))

        gt <- ggplot_gtable(ggplot_build(p))
        x_range = layer_scales(p)$x$range$range
        gt_widths=list(gt$widths[c(1:6)])

        # print(xy_list)
        # if (!is.na(vals$x)){
        #     p = p + geom_line(data=get_triangle(),
        #                       aes(x=x,y=y,group=group), inherit.aes=F)
        # }
        return(c(dt, dt_rev, center=center, max_color=max_color,
                 gt_widths=gt_widths))

    })


    updateVals <- function(x=NA, y=NA, start=NA, end=NA, region=NA){
        strand = ifelse(strand(vals$center)=='+', 1, -1)

        if (!is.na(region)){
            start = as.numeric(region[2]) - start(vals$center) * strand
            end = as.numeric(region[3]) - start(vals$center) * strand
        }
        if (is.na(x) & !is.na(start)){
            x = (start + end) / 2
            y = end - start
        } else{

            start = x - round(y/2)
            end = x + round(y/2)
        }
        if (is.na(region)){
            region = c(as.character(seqnames(vals$center)),
                       start(vals$center) + (start * strand),
                       start(vals$center) + (end * strand))
            if (strand==-1){
                region = region[c(1,3,2)]
            }
        }
        vals$x = x
        vals$y = y
        vals$start = start
        vals$end = end
        vals$fragment_warning = FALSE
        vals$region = region
    }

    get_triangle <- function(){
        triangle = data.table(x=c(vals$start, vals$x, vals$x, vals$end),
                              y=c(0, vals$y, vals$y, 0),
                              group= c(1,1,2,2))
        return(triangle)
    }

    observeEvent(input$plot_click, {
        print('click!')
        print(input$plot_click)
        x = round(input$plot_click$x)
        y = round(input$plot_click$y)
        updateVals(x=x, y=y)
        updateSelection()
        # input_list = dataInput()
        # triangle = get_triangle()
        # start = triangle[1,x]
        # end = triangle[4,x]
        # vals$score = input_list$flat_dt[x%in%start:end, sum(score)]
        # vals$score_rev = input_list$flat_dt_rev[x%in%start:end, sum(score)]
    })

    # observeEvent(input$plot_click_peak, {
    #     input_list = dataInput()
    #     nearPoints(input_list$peak_fwd, input$plot_click_peak)
    # })

    observeEvent(input$plot_brush, {

        ## weird translation necessary because grid.draw messes up coordinates
        ## had to print out all the coordinates of the grid lines
        ## and run a linear model.
        # dt = data.frame(x=c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1),
        #                 y=c(0.02180469, 0.1460047, 0.2783047, 0.3930547,
        #                     0.5159047, 0.6401047, 0.7643047, 0.8885047, 1.011355))
        # lm(dt)

        xmin = input$plot_brush$xmin * 1.01279 - 0.02432
        xmax = input$plot_brush$xmax * 1.01279 - 0.02432

        #
        # size = input$window[1] + input$window[2]
        #
        # xmin = xmin * size
        # xmax = xmax * size

        input_list = dataInput()

        x_size = -input$window[1] + input$window[2]

        x_start = xmin * x_size + input$window[1]
        x_end = xmax * x_size + input$window[1]

        updateVals(start=round(x_start), end=round(x_end))


        updateSelection()
    })


    output$frame <- renderPlotly({
        strand_color = c('+'='#228833', '-'='#EE6677')

        input_list = dataInput()
        if (!is.na(vals$center)){
            xregion <- getPlotRegion()
            tss_o = subjectHits(findOverlaps(xregion, tss_gr, ignore.strand=T))
            tss_info = as.data.table(tss_gr[tss_o])
            tss_info[, gene := factor(unlist(ensembl_dt[gsub('[.].*', '', gene_id),
                                                        'gene_name']))]
            tss_info[, sense := strand()]
            tss_info[, xstart := start(vals$center) - start]
            tss_info[, xend := ifelse(strand==strand(vals$center),
                                      xstart + 100, xstart-100)]
            tss_info[, y:=jitter(as.numeric(gene), amount=0.4)]

            tss_info[, y:=y[order(y, decreasing=T)],by=gene]

            tss_info[, color:=strand_color[strand]]
            tss_info[, xanchor:=ifelse(strand==strand(vals$center),
                                       'right', 'left')]

            tip = copy(tss_info)
            bottom = copy(tss_info)

            bottom[, y := 0]

            top = rbind(tip, tss_info)
            down = rbind(bottom, tss_info)

            tile_dt = data.table(gene=unique(tss_info$gene),
                                 xstart=input$window[1],
                                 xend=input$window[2],
                                 ystart=as.numeric(unique(tss_info$gene)) - 0.5,
                                 yend=as.numeric(unique(tss_info$gene)) + 0.5)
            grey_list = rep(c("grey65", "grey85"),10)

            p = ggplot(top, aes(y=y, x=xstart, colour = strand,
                                group=tx_name, label=gene, label2=gene_id,
                                label3=tx_name))  +
               geom_rect(data=tile_dt, inherit.aes=F,
                         aes(xmin=xstart, xmax=xend, ymin=ystart, ymax=yend,
                             fill=gene, alpha=0.3)) +
               # geom_text(data=tile_dt, inherit.aes=F,
               #           aes(x=xstart, y=as.numeric(label), label=label)) +
               scale_fill_manual(values=grey_list) +
               geom_line() +
               geom_line(data=down) +
               scale_color_manual(values=strand_color) +
               ggplot2::coord_cartesian(input$window, expand=0) +
               theme_bw() +
               xlab('Gencode v27 TSS annotations (liftOver to hg19)') +
               ggplot2::theme(axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank(),
                              axis.text.x=element_text(size=14))

            ply = ggplotly(p, tooltip=c('label', 'label2', 'label3')) %>%
               add_annotations(data=tss_info,
                               x = ~xend,
                               y = ~y,
                               xref = "x", yref = "y",
                               axref = "x", ayref = "y",
                               text = "",
                               showarrow = T,
                               xanchor=~xanchor,
                               arrowcolor= ~color,
                               ax = ~xstart,
                               ay = ~y) %>%
               layout(showlegend=F, margin = list(l=45, r=8)) %>%
               add_annotations(data=tile_dt, xanchor="left",xref='x',yref='y',
                               text=~gene, x=~xstart, y=~as.numeric(gene),
                               showarrow=F)
            ply
            # plot(JBrowseR("ViewHg19", location = loc))

        }
    })


    output$trianglePlot <- renderPlot({
        cutoff = ifelse(input$lib%in%c('HEPG2', 'K562'), input$cutoff42, input$cutoff23)
        start_time <- Sys.time()
        input_list = dataInput()

        end_time <- Sys.time()
        print(end_time - start_time)
        strand = strand(input_list$center)
        sec = start(input_list$center) * ifelse(strand=='+', 1, -1)

        p = plot_mat(input_list$mat_dt, input, cutoff, input_list$max_color,
                     vals$binsize) +
          ggtitle(paste("location on", seqnames(vals$center), '(hg19)')) +
          scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0))
        if (!is.na(vals$x)){
            p = p + geom_line(data=get_triangle(), aes(x=x,y=y,group=group),
                              inherit.aes=F)
        }
        gt <- ggplot_gtable(ggplot_build(p))
        p + ggplot2::theme(legend.position="none")

    })

    output$trianglePlot_rev <- renderPlot({
        input_list = dataInput()
        if (vals$show_minus==TRUE){
          cutoff = ifelse(input$lib%in%c('HEPG2', 'K562'), input$cutoff42, input$cutoff23)

          strand = strand(input_list$center)
          sec = start(input_list$center) * ifelse(strand=='+', 1, -1)

          p = plot_mat(input_list$mat_dt_rev, input, cutoff, input_list$max_color,
                       vals$binsize) +
            scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0))
          if (!is.na(vals$x)){
            p = p + geom_line(data=get_triangle(), aes(x=x,y=y,group=group),
                              inherit.aes=F)
          }
          gt <- ggplot_gtable(ggplot_build(p))
          p + ggplot2::theme(legend.position="none")
      }

    })


    output$legend <- renderPlot({
        cutoff = ifelse(input$lib%in%c('HEPG2', 'K562'), input$cutoff42, input$cutoff23)
        input_list = dataInput()
        p = plot_mat(input_list$mat_dt, input, cutoff, input_list$max_color) +
            labs(fill="predicted\nexpression")
        legend <- get_legend(p)
        grid.draw(legend)

    })

    output$flatPlot <- renderPlot({
        input_list = dataInput()
        strand = strand(input_list$center)
        sec = start(input_list$center) * ifelse(strand=='+', 1, -1)

        p2 = ggplot(input_list$flat_dt, aes(x=x_rel, y=score, color=score > 0)) +
          geom_histogram(stat='identity') +
          scale_fill_manual(values=c(T='blue',F='red')) +
          ggplot2::theme_bw() +
          xlab('position relative to point of interest') +
          ylab('coefficients') +
          ggplot2::theme(legend.position="none") +
          coord_cartesian(input$window) +
          scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0)) +
          ggplot2::theme(axis.text.y=element_text(size=14),
                         axis.title.x=element_text(size=14),
                         axis.title.y=element_text(size=14),
                         axis.text.x=element_text(size=14),
                         plot.margin=margin(t=0,r=20,b=0,l=0))
        if (!is.na(vals$x)){
          p2 = p2 + geom_vline(xintercept=vals$start) +
            geom_vline(xintercept=vals$end)
        }
        gt2 <- ggplot_gtable(ggplot_build(p2))
        gt2$widths[1:6] <- input_list$gt_widths
        grid.draw(gt2)
    })


    output$flatPlot_rev <- renderPlot({
        input_list = dataInput()
        if (vals$show_minus==TRUE){
            strand = strand(input_list$center)
            sec = start(input_list$center) * ifelse(strand=='+', 1, -1)
            seqname = seqnames(input_list$center)
            p2 = ggplot(input_list$flat_dt_rev, aes(x=x_rel, y=score, color=score > 0)) +
              geom_histogram(stat='identity') +
              scale_fill_manual(values=c(T='blue',F='red')) +
              ggplot2::theme_bw() +
              xlab('position relative to point of interest') +
              ylab('coefficients') +
              ggplot2::theme(legend.position="none") +
              coord_cartesian(input$window) +
              scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0)) +
              ggplot2::theme(axis.text.y=element_text(size=14),
                             axis.title.x=element_text(size=14),
                             axis.title.y=element_text(size=14),
                             axis.text.x=element_text(size=14),
                             plot.margin=margin(t=0,r=20,b=0,l=0))
            if (!is.na(vals$x)){
              p2 = p2 + geom_vline(xintercept=vals$start) +
                geom_vline(xintercept=vals$end)
            }
            gt2 <- ggplot_gtable(ggplot_build(p2))
            gt2$widths[1:6] <- input_list$gt_widths
            plot(gt2)
        }
    })


    output$text_minus <- renderText({
        if (vals$show_minus==TRUE){
            HTML("<h3>anti-sense orientation</h3>")
        }
    })

    output$text_plus <- renderText({
        if (!is.na(vals$center)){
            HTML("<h3>sense orientation</h3>")
        }
    })

    output$text_warning <- renderText({
        warning_vec = c()
        prefix = '<h3><p style="color:red;"><b>&lt;&lt;WARNING:'
        suffix = '&gt;&gt;</h3></b></p>'
        if (!is.na(vals$center) & vals$uptodate == FALSE){
            warning = paste(prefix, 'press submit to update configuration',
                            suffix)
            warning_vec = c(warning_vec, warning)
        }
        if (vals$fragment_warning == TRUE){
            warning = paste(prefix, 'please select a promoter fragment within',
                            'the region (or use the predict fragment tab)',
                            suffix)
            warning_vec = c(warning_vec, warning)
        }
        if (length(warning_vec) > 0){
            HTML(paste(warning_vec, sep='</br>'))
        }
    })

    # output$ROI_info <- renderText({
    #     steps = paste("<b><h3>How to:</h3>",
    #                   "<ol>",
    #                   "  <li>Enter a specific locus in the 'position of interest' field using:</li></b>",
    #                   "    <dd>- gene symbol (e.g. NUP214)",
    #                   "    <dd>- ensembl/gencode gene id (e.g. ENSG00000126883[.16])",
    #                   "    <dd>- ensemble/gencode transcript id (e.g. ENST00000359428[.5])",
    #                   "    <dd>- position on the genome in the form of [chromosome]:[position]:[strand] (used as center)</br></br>",
    #                   "  <b><li>Press Submit</li></br>")
    #
    #     steps = paste(steps,
    #                   "<li>Select a promoter fragment to predict its expression by:</li></b>",
    #                   "  <dd>- Clicking on a point within the triangle plot</dd>",
    #                   "  <dd>- Selecting a region in the coefficient plot <i>(click + drag)</i></dd>",
    #                   "  <dd>- Filling in promoter fragment coordinates <i>(and pressing update)</i></dd></br><b>")
    #
    #     HTML(paste0(steps, "</b>"))
    # })

    output$help_selection <- renderUI({
        if (!is.na(vals$center)){
            helpText(paste("expression prediction for selected promoter",
                           "fragment (for anti-sense use advanced options)"))
        }
    })

    output$selection <- renderTable({
        input_list = dataInput()

        if (!is.na(vals$x)){
            selection = vals$start:vals$end

            score = input_list$flat_dt[x_rel%in%selection, sum(score)]
            if (vals$show_minus){
                score_rev = input_list$flat_dt_rev[x_rel%in%selection, sum(score)]
                dt = data.table(prediction = exp(score),
                                prediction_rev = exp(score_rev))
                colnames(dt) = c('sense',
                                 'anti-sense')
            } else {
                dt = data.table(prediction = exp(score))
                colnames(dt) = c('sense')
            }
        } else {
            dt = data.table(prediction = NaN)
            colnames(dt) = c('sense')
        }
        return(dt)
    }, caption="Expression",
    caption.placement = getOption("xtable.caption.placement", "top"))

    parseLocation <- function(location){
        return(str_match(location, "(.*):([0-9]+)-([0-9]+)")[2:4])
    }

    getPlotRegion <- function(){
        xregion <- promoters(vals$center, upstream=input$window[1]*-1,
                             downstream=input$window[2])
        return(xregion)
    }

    observeEvent(input$update, {
        hg19_input = parseLocation(input$hg19_selection)
        hg38_input = parseLocation(input$hg38_selection)

        if (any(is.na(vals$region))){
            if (grepl('e.g.', hg19_input)[1]){
                hg19_input = liftOver(hg38_input)
            } else {
                hg38_input = liftOver(hg19_input, toHg19=F)
            }
        } else if (all(hg19_input==vals$region)){
            hg19_input = liftOver(hg38_input)
        } else {
            hg38_input = liftOver(hg19_input, toHg19=F)
        }

        xregion <- getPlotRegion()


        if (start(xregion) < hg19_input[2] & end(xregion) > hg19_input[3]){
            updateVals(region=hg19_input)
            updateSelection()

            input_list = dataInput()

        } else {
            vals$fragment_warning = TRUE
        }


        # vals$region = c(as.character(seqnames(center)),
        #                 start, end)
        #   d5 <<- editData(d5, info)
        # replaceData(proxy5, d5, resetPaging = FALSE)
    })

    output$ucsc <- renderText({
        input_list = dataInput()
        region = input_list$region

        location = paste0("position=", seqnames(region), "%3A", start(region), "-",
                          end(region))

        src = paste0(ucsc, paste(ucsc_opt, location, sep='&'))
        if (any(!is.na(vals$region))){
            jhighlight = paste0("highlight=", vals$region[1], "%3A", vals$region[2],
                                "-", vals$region[3])
            src = paste(src, jhighlight, sep='&')
        }
        HTML(paste0("<a class='btn btn-primary' href='", src,
                    "', target='_blank'>",
                    "open region in UCSC browser</a>"))
    })

    output$hg19_sel = renderUI({
        if (!is.na(vals$center)){
            fluidRow(textInput("hg19_selection", "promoter fragment hg19:",
                               value = "e.g. chr9:134000948-134001248"),
                     textInput("hg38_selection", "promoter fragment hg38:",
                               value = "e.g. chr9:134000948-134001248"),
                     actionButton("update", "Update"))

        }
    })

    fragmentInput <- observeEvent(input$fragment_file,{
        upload = paste(readLines(input$fragment_file$datapath), collapse="\n")
        updateTextAreaInput(session, "text_fragments", value = upload)
    })

    fragmentInput <- eventReactive(input$go_fragments, {
        width_offset = offset_coefs$width[[1]]

        frag_dt = fread(input$text_fragments, stringsAsFactors=F,
                        fill=T, header=F, sep='\t',
                        sep2=' ')
        colnames(frag_dt) = c('seqnames', 'start', 'end',
                              'name', 'color', 'strand')[1:ncol(frag_dt)]
        frag_dt[,strand:=ifelse(strand%in%c('+','-'), strand, "*")]
        hg19_dt = copy(frag_dt)
        if (input$hg_version=='hg38'){
           suppressWarnings(hg19_dt[, c('seqnames', 'start', 'end') :=
                                      as.list(liftOver(c(seqnames, start, end))),
                                    by=eval(colnames(frag_dt))])
        }

        frag_gr = makeGRangesFromDataFrame(hg19_dt)

        track_names = c('K562_plus', 'K562_minus', 'HEPG2_plus',
                        'HEPG2_minus')

        score_list = lapply(track_names, function(track){
            file = paste0('data/',
                          'SURE_elasticNet_allele_',
                          track, '.bw')
            glm <- import(file, selection = BigWigSelection(frag_gr),
                          as="NumericList")
            score = lapply(glm, function(x){exp(sum(x + width_offset))})
            return(unlist(score))
        })
        score_dt = as.data.table(do.call(cbind, score_list))
        colnames(score_dt) = track_names

        use_strand <- function(score_dt, strand){
          if (strand=="+"){
            dt = score_dt
            colnames(dt) = gsub("_plus", "_sense", colnames(dt))
            colnames(dt) = gsub("_minus", "_antisense", colnames(dt))
          } else if (strand=="-"){
            dt = score_dt
            colnames(dt) = gsub("_minus", "_sense", colnames(dt))
            colnames(dt) = gsub("_plus", "_antisense", colnames(dt))
            i_vec = unlist(lapply(1:(ncol(dt)/2), function(x){x* 2:1}))
            dt = dt[,i_vec,with=F]
          } else {
            dt = score_dt
          }
          return(dt)
        }

        score_dt[,strand:=as.character(strand(frag_gr))]
        if ("strand" %in% colnames(frag_dt)){
          score_dt = score_dt[,use_strand(.SD, strand), by="strand"]
        }
        result = cbind(frag_dt, score_dt[,-"strand"])

        return(result)
    })

    output$download <- renderUI({
        if (input$go_fragments!=0){
          fluidRow(
            downloadButton("download_csv", "Download as csv"),
            downloadButton("download_tsv", "Download as tsv")
          )
        }
    })

    output$fragment_result <- renderDataTable({

        return(fragmentInput())
    })

    output$download_csv <- downloadHandler(
      filename = function() {
        paste("SuRE-GLM_predictions-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        fwrite(fragmentInput(), file, sep=',')
      }
    )

    output$download_tsv <- downloadHandler(
      filename = function() {
        paste("SuRE-GLM_predictions-", Sys.Date(), ".tsv", sep="")
      },
      content = function(file) {
        fwrite(fragmentInput(), file, sep='\t')
      }
    )
})



# shinyApp(ui = ui, server = server, options = list(display.mode="showcase"))
