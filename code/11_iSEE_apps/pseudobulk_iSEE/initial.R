initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "PCA", XAxis = 1L, YAxis = 2L,
                                          FacetRowByColData = "brnum", FacetColumnByColData = "brnum",
                                          ColorByColumnData = "domain", ColorByFeatureNameAssay = "logcounts",
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "brnum",
                                          SizeByColumnData = "age", TooltipColumnData = character(0),
                                          FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data",
                                          ColorByDefaultColor = "#000000", ColorByFeatureName = "LINC01409",
                                          ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                          ColorBySampleName = "1", ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
                                          ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
                                          ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = TRUE,
                                          VisualChoices = c("Color", "Size"), ContourAdd = FALSE, ContourColor = "#0000FF",
                                          PointSize = 2, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
                                          CustomLabels = FALSE, CustomLabelsText = "1", FontSize = 1,
                                          LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
                                          LabelCenters = FALSE, LabelCentersBy = "brnum", LabelCentersColor = "#000000",
                                          VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                              "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 4L,
                                          SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
                                          DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                          RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                          SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE,
                                        CustomRowsText = "PPFIA2\nAMPH\nFNDC1\nGFRA1\nKRT17\nMEF2C\nGAD2\nAPOC1\nMT1G\nSLC1A2\nSLC1A3\nSFRP2\nMOBP\nACTA2\nPRLR",
                                        ClusterRows = FALSE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = c("domain",
                                                                                                           "broad.domain"), RowData = character(0), CustomBounds = FALSE,
                                        LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = TRUE,
                                        AssayScaleRows = TRUE, DivergentColormap = "blue < white < red",
                                        ShowDimNames = "Rows", LegendPosition = "Bottom", LegendDirection = "Horizontal",
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10,
                                        ShowColumnSelection = FALSE, OrderColumnSelection = TRUE,
                                        VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                            "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 4L,
                                        SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
                                        RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                        RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                        SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
                                      XAxisColumnData = "domain", XAxisFeatureName = "LINC01409",
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
                                      YAxisFeatureName = "CLSTN3", YAxisFeatureSource = "---",
                                      YAxisFeatureDynamicSource = FALSE, FacetRowByColData = "brnum",
                                      FacetColumnByColData = "brnum", ColorByColumnData = "domain",
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
                                      ShapeByColumnData = "brnum", SizeByColumnData = "age", TooltipColumnData = character(0),
                                      FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data",
                                      ColorByDefaultColor = "#000000", ColorByFeatureName = "LINC01409",
                                      ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                      ColorBySampleName = "1", ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
                                      ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
                                      ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = TRUE,
                                      VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
                                      PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
                                      CustomLabels = FALSE, CustomLabelsText = "1", FontSize = 1,
                                      LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
                                      LabelCenters = FALSE, LabelCentersBy = "brnum", LabelCentersColor = "#000000",
                                      VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                          "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 4L,
                                      SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
                                      DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                      RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                      SelectionHistory = list())
