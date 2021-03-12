library(rtracklayer)
library(data.table)
library(shiny)
library(AnnotationDbi)
library(ggplot2)
library(cowplot)
library(grid)
library(stringr)
library(here)


peak_file  = c('K562' = 'data/SuRE_elasticNet_peaks_K562.bed',
               'HEPG2' = 'data/SuRE_elasticNet_peaks_K562.bed')

rd_file = 'data/promoter_triangle_input.RData'

offset_coefs = readRDS("data/optimalFitOffsetCoefficients_K562_spatial.rds")

jbrowse = "***REMOVED***jbrowse-hg19/index.html?tracklist=0&nav=0&overview=0"
jtrack_start = paste("tracks=gencode.v27lift37",
                     "GSM1480321_K562_GROcap_wTAP_minus",
                     "GSM1480321_K562_GROcap_wTAP_plus",
                     "TTseq_K562_rep2_plus_hg19",
                     "TTseq_K562_rep2_minus_hg19", sep='%2C')

hg38ToHg19 = "data/hg38ToHg19.over.chain"
hg19ToHg38 = "data/hg19ToHg38.over.chain"


load(rd_file)

get_center <- function(text_input, lookup_list, tss_gr){
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
  } else if (grepl('ENST[0-9]+[.][0-9]+', text_input)){
    center = tss_gr[text_input]
  } else {
    i = 1
    is_match = F
    while (!is_match==T) {
      is_match = grepl(match_list[i], text_input)
      match_name = names(match_list)[i]
      i = i + 1
    }
    transcript_id = lookup_list[[match_name]][text_input, transcript_id]
    center = tss_gr[transcript_id]
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
    y_mat <- abs(end_mat - start_mat)[!is.na(mat)]
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
    if (strand=='plus'){
        dt = data.table(x=x_mat, y=y_mat, score=mat[!is.na(mat)],
                      x_rel = (x_mat - center_pos))
        mat_dt = dt[x_rel <= downstream & x_rel >= -upstream, ]

        flat_dt = data.table(score=glm,
                             x = start(region):end(region))
        flat_dt[,x_rel := x - center_pos]
    } else {
        dt = data.table(x=x_mat, y=y_mat, score=mat[!is.na(mat)],
                      x_rel = (center_pos - x_mat))
        mat_dt = dt[x_rel <= upstream & x_rel >= -downstream, ]


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

colorlut <- colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
                               "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61",
                               "#F46D43", "#D73027", "#B91326",
                               "#A50026"))(101) # height color lookup table

plot_mat <-function(mat, input, cutoff, max_color){
    print(exp(max_color))
    print(max(exp(mat$score)))
  p = ggplot(mat, aes(x=x_rel, y=y, fill=exp(score))) +
    geom_tile(width=input$binsize, height=input$binsize) +
    theme_bw() +
    scale_fill_gradientn(colours=colorlut, limits=c(0,exp(max_color))) +
    coord_cartesian(ylim=c(0,cutoff)) +
    ylab('reporter length') +
    xlab('##chromosome##') +
    coord_fixed() +
    theme(axis.title.x=element_blank())
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
    chain = ifelse (toHg19, hg38ToHg19, hg19ToHg38)
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
        center = NA
    )

    dataInput <- eventReactive(input$go, {
        center = get_center(input$ROI, lookup_list, tss_gr)
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
                         cutoff = cutoff, binsize=input$binsize, offset=offset,
                         filePart = lib_vec[input$lib])

        dt_rev = triangle_dt(center, upstream=input$window[1]*-1,
                             downstream=input$window[2],
                             cutoff = cutoff, binsize=input$binsize, offset=offset,
                             filePart = lib_vec[input$lib], rev_comp = T)
        names(dt_rev) = paste0(names(dt_rev), '_rev')
        strand = strand(center)
        sec = start(center) * ifelse(strand=='+', 1, -1)


        peak_dt = suppressWarnings(get_peaks(dt$region,
                                             peak_file[input$lib]))
        if (nrow(peak_dt)==0){
            peak_melt = peak_dt
        } else {
            peak_dt[,group:=1:nrow(peak_dt)]
            peak_melt = melt(peak_dt, measure.vars=c('start', 'end'),
                             value.name='x')
            if (as.character(strand) == '+'){
                peak_melt[,x_rel := x - start(center)]
                peak_dt[, strand_order:= ifelse(strand=='+',2,1)]
                max_peak = peak_dt[order(strand_order,score, decreasing=T), ][1,]
                updateVals(start=max_peak[['start']] - start(center),
                           end=max_peak[['end']] - start(center))
            } else {
                peak_melt[,x_rel := start(center) - x]
                peak_dt[, strand_order:= ifelse(strand=='+',1,2)]
                max_peak = peak_dt[order(strand_order,score, decreasing=T), ][1,]

                updateVals(start=start(center) - max_peak[['start']],
                           end=start(center) - max_peak[['end']])
            }
        }
        max_color = max(dt$mat_dt[, max(score)],
                        dt_rev$mat_dt_rev[, max(score)])

        p = plot_mat(dt$mat_dt, input, cutoff, max_color) +
          scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0))
        # print(xy_list)
        # if (!is.na(vals$x)){
        #     p = p + geom_line(data=get_triangle(),
        #                       aes(x=x,y=y,group=group), inherit.aes=F)
        # }
        gt <- ggplot_gtable(ggplot_build(p))
        x_range = layer_scales(p)$x$range$range
        return(c(dt, dt_rev, center=center, max_color=max_color,
                 gt_widths=list(gt$widths[c(1:6)]),
                 list(peak_fwd=peak_melt[strand==as.character(strand(center)), ],
                      peak_rev=peak_melt[strand!=as.character(strand(center)), ])))
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
            print(x)
            print(y)
            start = x - round(y/2)
            end = x + round(y/2)
        }
        if (is.na(region)){
            region = c(as.character(seqnames(vals$center)),
                       start(vals$center) + (start * strand),
                       start(vals$center) + (end * strand))
        }
        vals$x = x
        vals$y = y
        vals$start = start
        vals$end = end
        print(vals)
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
        xmin = round(input$plot_brush$xmin)
        xmax = round(input$plot_brush$xmax)

        left = input$plot_brush$domain$left
        right = input$plot_brush$domain$right

        xmin = xmin * (right - left)
        xmax = xmax * (right - left)

        print(input$plot_brush)
        input_list = dataInput()

        x_size = input$window[2] - input$window[1]
        x_start = input$window[1] + xmin * x_size
        x_end = input$window[1] + xmax * x_size

        updateVals(start=x_start, end=x_end)
        #
        # y = x_end - x_start
        # x = (x_start + x_end) / 2
        # vals$triangle <- data.frame(x=c(x_start, x, x, x_end),
        #                             y=c(0, y, y, 0),
        #                             group= c(1,1,2,2))
        #
        # input_list$x_range[1]
        # center = var$center
        # start = start(center) + x_start
        # end = start(center) + x_end

        updateSelection()
        # vals$score = input_list$flat_dt[x%in%start:end, sum(score)]
        # vals$score_rev = input_list$flat_dt_rev[x%in%start:end, sum(score)]
    })


    output$trianglePlot <- renderPlot({

        cutoff = ifelse(input$lib%in%c('HEPG2', 'K562'), input$cutoff42, input$cutoff23)
        start_time <- Sys.time()
        input_list = dataInput()

        end_time <- Sys.time()
        print(end_time - start_time)
        strand = strand(input_list$center)
        sec = start(input_list$center) * ifelse(strand=='+', 1, -1)

        p = plot_mat(input_list$mat_dt, input, cutoff, input_list$max_color) +
          scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0))
        if (!is.na(vals$x)){
            p = p + geom_line(data=get_triangle(), aes(x=x,y=y,group=group),
                              inherit.aes=F)
        }
        gt <- ggplot_gtable(ggplot_build(p))
        plot(p + theme(legend.position="none"))
    })

    output$trianglePlot_rev <- renderPlot({

      cutoff = ifelse(input$lib%in%c('HEPG2', 'K562'), input$cutoff42, input$cutoff23)

      input_list = dataInput()
      strand = strand(input_list$center)
      sec = start(input_list$center) * ifelse(strand=='+', 1, -1)

      p = plot_mat(input_list$mat_dt_rev, input, cutoff, input_list$max_color) +
        scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0))
      if (!is.na(vals$x)){
        p = p + geom_line(data=get_triangle(), aes(x=x,y=y,group=group),
                          inherit.aes=F)
      }
      gt <- ggplot_gtable(ggplot_build(p))
      p + theme(legend.position="none")
    })


    output$legend <- renderPlot({

        cutoff = ifelse(input$lib%in%c('HEPG2', 'K562'), input$cutoff42, input$cutoff23)
        input_list = dataInput()
        p = plot_mat(input_list$mat_dt, input, cutoff, input_list$max_color) +
            theme(legend.title=element_blank())
        legend <- get_legend(p)
        plot(legend)
    })

    output$flatPlot <- renderPlot({
        input_list = dataInput()
        strand = strand(input_list$center)
        sec = start(input_list$center) * ifelse(strand=='+', 1, -1)

        p2 = ggplot(input_list$flat_dt, aes(x=x_rel, y=score, color=score > 0)) +
          geom_histogram(stat='identity') +
          scale_fill_manual(values=c(T='blue',F='red')) +
          theme_bw() +
          xlab('coefficients') +
          theme(legend.position="none") +
          coord_cartesian(input$window) +
          scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0))
        if (!is.na(vals$x)){
          p2 = p2 + geom_vline(xintercept=vals$start) +
            geom_vline(xintercept=vals$end)
        }
        gt2 <- ggplot_gtable(ggplot_build(p2))
        gt2$widths[1:6] <- input_list$gt_widths
        grid.draw(gt2)
    })
    #
    # output$peakPlot <- renderPlot({
    #     input_list = dataInput()
    #     strand = strand(input_list$center)
    #     sec = start(input_list$center) * ifelse(strand=='+', 1, -1)
    #     peak_dt = input_list$peak_fwd
    #     print(peak_dt)
    #     if (nrow(peak_dt) > 0){
    #         p2 = ggplot(peak_dt, aes(x=x_rel, y=factor(group),
    #                                  color=exp(score))) +
    #           geom_line(size=3) +
    #           scale_color_gradientn(colours=colorlut,
    #                                 limits=c(0, exp(input_list$max_color))) +
    #           theme_bw() +
    #           xlab('peaks') +
    #           coord_cartesian(input$window) +
    #           scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0)) +
    #           scale_y_discrete(breaks=NULL) +
    #           theme(axis.title.y=element_blank(),
    #                 legend.position="none")
    #         if (!is.na(vals$x)){
    #           p2 = p2 + geom_vline(xintercept=vals$start) +
    #             geom_vline(xintercept=vals$end)
    #         }
    #         gt2 <- ggplot_gtable(ggplot_build(p2))
    #         gt2$widths[1:6] <- input_list$gt_widths
    #     } else {
    #         gt2 = 0
    #     }
    #     plot(gt2)
    # })

    output$flatPlot_rev <- renderPlot({
        input_list = dataInput()
        strand = strand(input_list$center)
        sec = start(input_list$center) * ifelse(strand=='+', 1, -1)

        p2 = ggplot(input_list$flat_dt_rev, aes(x=x_rel, y=score, color=score > 0)) +
          geom_histogram(stat='identity') +
          scale_fill_manual(values=c(T='blue',F='red')) +
          theme_bw() +
          xlab('coefficients') +
          theme(legend.position="none") +
          coord_cartesian(input$window) +
          scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0))
        if (!is.na(vals$x)){
          p2 = p2 + geom_vline(xintercept=vals$start) +
            geom_vline(xintercept=vals$end)
        }
        gt2 <- ggplot_gtable(ggplot_build(p2))
        gt2$widths[1:6] <- input_list$gt_widths
        plot(gt2)
    })
    #
    # output$peakPlot_rev <- renderPlot({
    #     input_list = dataInput()
    #     strand = strand(input_list$center)
    #     sec = start(input_list$center) * ifelse(strand=='+', 1, -1)
    #     peak_dt = input_list$peak_rev
    #     print(peak_dt$score)
    #     print(peak_dt[,exp(score)])
    #     if (nrow(peak_dt) > 0){
    #         p2 = ggplot(peak_dt, aes(x=x_rel, y=factor(group),
    #                                  color=exp(score))) +
    #           geom_line(size=3) +
    #           scale_color_gradientn(colours=colorlut,
    #                                 limits=c(0, exp(input_list$max_color))) +
    #           theme_bw() +
    #           xlab('peaks') +
    #           coord_cartesian(input$window) +
    #           scale_x_continuous(sec.axis=~abs(. + sec), expand=c(0,0)) +
    #           scale_y_discrete(breaks=NULL) +
    #           theme(axis.title.y=element_blank(),
    #                 legend.position="none")
    #         if (!is.na(vals$x)){
    #           p2 = p2 + geom_vline(xintercept=vals$start) +
    #             geom_vline(xintercept=vals$end)
    #         }
    #         gt2 <- ggplot_gtable(ggplot_build(p2))
    #         gt2$widths[1:6] <- input_list$gt_widths
    #     } else {
    #         gt2 = 0
    #     }
    #     plot(gt2)
    # })


    output$selection <- renderTable({
        input_list = dataInput()

        if (!is.na(vals$x)){
            selection = vals$start:vals$end

            score = input_list$flat_dt[x_rel%in%selection, sum(score)]
            score_rev = input_list$flat_dt_rev[x_rel%in%selection, sum(score)]
            dt = data.table(prediction = exp(score),
                            prediction_rev = exp(score_rev))
            colnames(dt) = c('sense expression',
                             'anti-sense expression')
        } else {
            dt = data.table(prediction = NaN,
                            prediction_rev = NaN)
            colnames(dt) = c('sense expression',
                             'anti-sense expression')
        }
        return(dt)
    })

    parseLocation <- function(location){
        return(str_match(location, "(.*):([0-9]+)-([0-9]+)")[2:4])
    }

    observeEvent(input$update, {
        hg19_input = parseLocation(input$hg19_selection)
        hg38_input = parseLocation(input$hg38_selection)

        if (all(hg19_input==vals$region)){
            hg19_input = liftOver(hg38_input)
        } else {
            hg38_input = liftOver(hg19_input, toHg19=F)
        }
        updateVals(region=hg19_input)
        updateSelection()
        input_list = dataInput()


        # vals$region = c(as.character(seqnames(center)),
        #                 start, end)
        #   d5 <<- editData(d5, info)
        # replaceData(proxy5, d5, resetPaging = FALSE)
    })

    output$jbrowse <- renderUI({
        input_list = dataInput()
        region = input_list$region

        jloc = paste0("loc=", seqnames(region), "%3A", start(region), "...",
                      end(region))
        sure_tracks = paste0(lib_vec[input$lib], c('plus', 'minus'))
        jtracks = paste(c(jtrack_start, sure_tracks), collapse="%2C")
        src = paste(jbrowse, jloc, jtracks, sep='&')
        if (any(!is.na(vals$region))){
            jhighlight = paste0("highlight=", vals$region[1], "%3A", vals$region[2],
                                "...", vals$region[3])
            src = paste(src, jhighlight, sep='&')
        }
        my_test <- tags$iframe(src=src, height=200, width=800,
                               onload="vertical_scrollbar.scrollTo(0,0);")
        my_test
    })
})
# shinyApp(ui = ui, server = server, options = list(display.mode="showcase"))
