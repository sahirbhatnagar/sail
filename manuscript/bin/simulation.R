## ---- simulation-results ----

df <- readRDS("/home/sahir/git_repositories/sail/my_sims/simulation_results/apr_25_2018_results.rds")
df <- df %>% separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
                      sep = "/")

DT <- as.data.table(df, stringsAsFactors = FALSE)

# DT[parameterIndex=="parameterIndex_1", table(Method)] %>% names %>% dput
DT[Method=="Adaptivesail", Method := "Asail"]

appender <- function(string) TeX(paste(string))

DT[, method:=factor(Method, levels = c("lasso", "lassoBT", "GLinternet", "HierBasis", "SPAM", "gamsel",
                                       "sail", "Asail"))]
# DT[, table(method)]
DT[, scenario:= as.numeric(as.character(stringr::str_extract_all(parameterIndex, "\\d", simplify = T)))]
DT[, scen:=ifelse(scenario==1,"Strong Hierarchy",ifelse(scenario==2, "Weak Hierarchy", ifelse(scenario==3,"Interactions Only",ifelse(scenario==4, "Strong Hierarchy (Linear)", "Main Effects Only"))))]
DT[, scen:=factor(scen, levels = c("Strong Hierarchy", "Weak Hierarchy","Interactions Only","Strong Hierarchy (Linear)", "Main Effects Only"))]
# DT$scen %>% table
#Truth obeys strong hierarchy (parameterIndex = 1)
#Truth obeys weak hierarchy (parameterIndex = 2)
#Truth only has interactions (parameterIndex = 3)
#Truth is linear (parameterIndex = 4)
#Truth only has main effects (parameterIndex = 5)


p1_mse <- ggplot(DT, aes(method, mse, fill = method)) +
    ggplot2::geom_boxplot() +
    # gg_sy +
    # facet_rep_wrap(~scen, scales = "free", ncol = 2,
    #                repeat.tick.labels = 'left',
    #                labeller = as_labeller(appender,
    #                                       default = label_parsed)) +
    facet_rep_wrap(~scen, scales = "free", ncol = 2,
                   repeat.tick.labels = 'left',
                   labeller = as_labeller(appender,
                                          default = label_parsed)) +
    scale_fill_manual(values=cbbPalette, guide=guide_legend(ncol=2)) +
    # ggplot2::labs(y = "Test Set MSE", title = "") + xlab("") +
    labs(x="", y="Test Set MSE",
         title="Simulation Study Results: Test Set MSE",
         subtitle="Based on 200 simulations",
         caption="") +
    # panel_border()+
    # background_grid()+
    theme_ipsum_rc() + theme(legend.position = "right", axis.text.x = element_text(angle = 20, hjust = 1))

# , legend.text=element_text(size=18)


reposition_legend(p1_mse, 'center', panel='panel-2-3')


