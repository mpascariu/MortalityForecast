# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# Last update: Tue Mar 26 16:49:05 2019
# --------------------------------------------------- #
remove(list = ls())
library(xtable)
library(MortalityForecast)
library(magrittr)
load("dxForecastResults.Rdata")
options(xtable.comment = FALSE)

mixRows <- function(x, y) paste0(format(round(x, 2), nsmall = 2), " (", y, ")")

summary_table <- function(A, R, ns, period, label = "") {
  A$Model <- change_model_factor_levels(A$Model)
  R$Model <- change_model_factor_levels(R$Model)

  AR <- A
  for (i in 4:9) {
    a <- unlist(A[, i])
    r <- unlist(R[, i])
    AR[, i] <- mixRows(a, r)

  }

  AR$GC <- paste0("(", R$ModelRanking, ")")

  # Summary table
  if (missing(ns)) ns <- "multiple"
  if (missing(period)) period <- "observed"

  print(xtable(AR[AR$Scenario == "Total", c("Model", "ME", "MAE", "MAPE",
                                            "sMAPE", "sMRAE", "MASE", "GC")],
               caption = paste0("Forecast accuracy measures aggregated over ",
                                ns, " scenarios in the ", period, " period"),
               label = paste0("tbl:res_", label),
               align = "llccccccc"),
        table.placement = "!ht",
        caption.placement = "top",
        scalebox = 0.9,
        include.rownames = FALSE,
        booktabs = TRUE)
}


tex_intro1 <- c("\\documentclass[T0_MEM]{subfiles}\n",
               "\\linespread{1.1}\n",
               "\\begin{document}\n\n",
               "\\section*{Tables}\n\n")
tbl_EW <-  summary_table(A[["GBRTENW"]], R[["GBRTENW"]], ns = 18,
                         "1960--2016", label = "GBRTENW")  # EW
tbl_MEM <- summary_table(A2, R2, ns = 18, "1960--2016", label = "MEM") # EW - MEM only
tex_end1 <- c("\n\\end{document}")


# ---------------------------------------------------------

M <- as.data.frame(matrix(NA, nrow = 10, ncol = 6,
                              dimnames = list(cntrN,
                                              c("MRWD",
                                                "Lee-- Carter",
                                                "Hyndman-- Ullah",
                                                "Renshaw-- Haberman",
                                                "Oepppen",
                                                "MEM"))))

for (i in 1:10) M[i, ] <- unlist(R[[i]][1:6, "ModelRanking"]) %>% paste0("(", ., ")")


tbl_ranks <- print(xtable(M,
                          caption = paste("General classification of the model",
                                          "predictive power for male mortality",
                                          "experience in England and Wales and",
                                          "nine developed countries.",
                                          "Period analised: 1960--2016."),
                          label = "tbl:rank_countries",
                          align = "lp{1.5cm}p{1.5cm}p{2.5cm}p{2.5cm}p{1.5cm}p{1.5cm}"),
                   table.placement = "!ht",
                   caption.placement = "top",
                   scalebox = 0.9,
                   include.rownames = TRUE,
                   booktabs = TRUE)


Conn1 <- file("T3_Tables.tex")
writeLines(c(tex_intro1, tbl_EW, tbl_MEM, tbl_ranks, tex_end1),
           Conn1, sep = "")
close(Conn1)

# ----------------------------------------------
# Generate Appendix Tex file

C_Appendix <- cntr[cntr != "GBRTENW"]
CN_Appendix <- cntrN[cntr != "GBRTENW"]



tex_intro <- c("\\documentclass[T0_MEM]{subfiles}\n",
               "\\linespread{1.1}\n",
               "\\begin{document}\n\n",
               "\\section{Appendix B -- Out-of-sample forecasts}\n",
               "\\label{sec:more_countries}\n",
               "In this section we present the results of our study
               performed in nine additional developed countries.\n\n")
tex_body <- c()
tex_end <- c("\\newpage\n",
             "\\end{document}")

for (i in seq_along(C_Appendix)) {
  K  <- C_Appendix[i]
  KN <- CN_Appendix[i]
  ns <- as.numeric(rev(A[[K]]$Scenario)[i])

  Bi <- B[[K]]
  ns <- nrow(Bi$scenarios) # no of scenarios
  last_scenario <- Bi$results[[ns]]
  y_fit <- last_scenario$input$y.fit
  y_for <- last_scenario$input$y.for
  scenario <- paste(min(y_fit), max(y_fit), max(y_for), sep = "--")
  y_range <- paste(range(Bi$input$y), collapse = "--")


  tex0 <-  c("\\newpage % ---------------------------\n",
             "\\subsection{Out-of-sample forecasts: ",
             KN, ", Male population, ", y_range, "}","\n")
  texTable <- summary_table(A[[K]], R[[K]], ns, period = y_range, label = K)
  texFigures <- c(
    "\n",
    "\\begin{figure}[!hb]", "\n",
    "  \\includegraphics[width=1\\linewidth]{figure/Figure_", K,"_ex}",  "\n",
    "  \\includegraphics[width=1\\linewidth]{figure/Figure_", K, "_mx}",  "\n",
    "  \\caption{Out-of-sample forecast of the remaining life expectancy ",
    "and central death rates at various ages using the five mortality models (",
     KN,", Male population, Scenario ", ns, ": ", scenario, ")}", "\n",
    "\\end{figure}",
    "\n\n"
  )
  tex_body <- c(tex_body, tex0, texTable, texFigures)
}

tex_full <- c(tex_intro, tex_body, tex_end)

fileConn <- file("T4_Appendix_B.tex")
writeLines(tex_full, fileConn, sep = "")
close(fileConn)



