initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "UMAP", XAxis = 1L, YAxis = 2L,
                                          FacetRowByColData = "Sample", FacetColumnByColData = "Sample",
                                          ColorByColumnData = "fine.cell.class", ColorByFeatureNameAssay = "logcounts",
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample",
                                          SizeByColumnData = "sum", TooltipColumnData = character(0),
                                          FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data",
                                          ColorByDefaultColor = "#000000", ColorByFeatureName = "TTR",
                                          ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                          ColorBySampleName = "1_AAACCCACAACGATCT-1", ColorBySampleSource = "---",
                                          ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                          SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
                                          VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
                                          ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
                                          Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
                                          CustomLabelsText = "1_AAACCCACAACGATCT-1", FontSize = 1,
                                          LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
                                          LabelCenters = FALSE, LabelCentersBy = "Sample", LabelCentersColor = "#000000",
                                          VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                              "numeric_version"))), PanelId = 1L, PanelHeight = 500L, PanelWidth = 4L,
                                          SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
                                          DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                          RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                          SelectionHistory = list())


################################################################################
# Settings for Complex heatmap 2
################################################################################

initial[["ComplexHeatmapPlot2"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE,
                                        CustomRowsText = "CLSTN3\nC3\nGFAP\nMOBP\nPDGFRA\nDNAH11\nEBF1",
                                        ClusterRows = FALSE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = "fine.cell.class",
                                        RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_,
                                        UpperBound = NA_real_, AssayCenterRows = FALSE, AssayScaleRows = FALSE,
                                        DivergentColormap = "purple < black < yellow", ShowDimNames = "Rows",
                                        LegendPosition = "Right", LegendDirection = "Horizontal",
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10,
                                        ShowColumnSelection = FALSE, OrderColumnSelection = TRUE,
                                        VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                            "numeric_version"))), PanelId = 2L, PanelHeight = 500L, PanelWidth = 6L,
                                        SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
                                        RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                        RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                        SelectionHistory = list())
