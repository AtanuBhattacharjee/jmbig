# ============================================================================
# jmBIG dynamic prediction — deployable Shiny app (single file, self-contained)
# Simulate a joint-model dataset, fit jmBIG::jmbayesBig, then show individual
# dynamic prediction (biomarker forecast + conditional survival) with 95%
# credible bands and the patient's forecast error.
# ============================================================================
library(shiny); library(jmBIG); library(JMbayes2); library(nlme); library(survival)

sim <- function(n, seed = 1, lambda = 0.1, alpha = 0.5, tmax = 8) {
  set.seed(seed); D <- matrix(c(1, .1, .1, .3), 2); x1 <- rbinom(n, 1, .5)
  b <- MASS::mvrnorm(n, c(0, 0), D); e0 <- 5 + b[, 1] + .5 * x1; e1 <- -.3 + b[, 2]
  A <- lambda * exp(.3 * x1 + alpha * (e0 - 5)); bc <- alpha * e1; u <- runif(n)
  Te <- numeric(n); sm <- abs(bc) < 1e-8; Te[sm] <- -log(u[sm]) / A[sm]
  ar <- 1 + (-log(u[!sm]) * bc[!sm]) / A[!sm]; Te[!sm] <- ifelse(ar > 0, log(ar) / bc[!sm], Inf)
  ot <- pmin(Te, tmax); st <- as.integer(Te <= tmax); g <- seq(0, tmax, .5)
  L <- do.call(rbind, lapply(1:n, function(i) {
    v <- g[g <= ot[i]]; if (!length(v)) v <- 0
    data.frame(id = i, visit = v, x1 = x1[i], y = e0[i] + e1[i] * v + rnorm(length(v), 0, .6))
  }))
  list(dl = L, ds = data.frame(id = 1:n, time = ot, status = st, x1 = x1), tmax = tmax)
}

fitbig <- function(dl, ds)
  jmbayesBig(data.frame(dl), data.frame(ds), y ~ visit + x1,
             Surv(time, status) ~ x1, rd = ~ visit | id, timeVar = "visit",
             id = "id", samplesize = 200, nchain = 1, niter = 1000, nburnin = 400)

nd_of <- function(dl, sid, upto = Inf) {
  nd <- subset(dl, id == sid & visit <= upto); nd$time <- max(nd$visit); nd$status <- 0; nd
}
patn <- function(o, u) {
  d <- as.data.frame(o); p <- grep("^pred_", names(d), value = TRUE)[1]
  oo <- order(d$visit); approx(d$visit[oo], d[[p]][oo], u, rule = 2)$y
}

ui <- fluidPage(
  titlePanel("jmBIG — individual dynamic prediction"),
  sidebarLayout(
    sidebarPanel(
      helpText("1. Simulate a dataset, 2. Fit the model, 3. Explore a patient."),
      sliderInput("n", "patients", 100, 1000, 300, 100),
      sliderInput("alpha", "association (alpha)", -1, 1, 0.5, 0.1),
      actionButton("gen", "Generate data", class = "btn-primary"),
      tags$hr(),
      actionButton("fit", "Fit model", class = "btn-success"),
      tags$hr(),
      uiOutput("pid"), uiOutput("lm")
    ),
    mainPanel(
      plotOutput("pl", height = "360px"),
      tags$h4("Individual prediction error"),
      verbatimTextOutput("err")
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(dl = NULL, ds = NULL, big = NULL, tmax = 8)

  observeEvent(input$gen, {
    d <- sim(input$n, alpha = input$alpha); rv$dl <- d$dl; rv$ds <- d$ds
    rv$tmax <- d$tmax; rv$big <- NULL
    showNotification(sprintf("n=%d, events=%.0f%%", input$n, 100 * mean(d$ds$status)))
  })
  observeEvent(input$fit, {
    req(rv$dl); withProgress(message = "Fitting jmBIG…", value = .5, {
      rv$big <- fitbig(rv$dl, rv$ds) }); showNotification("fitted")
  })
  output$pid <- renderUI({ req(rv$ds); selectInput("pid", "patient id", sort(unique(rv$ds$id))) })
  output$lm <- renderUI({
    req(rv$dl, input$pid)
    v <- sort(unique(subset(rv$dl, id == as.numeric(input$pid))$visit))
    if (length(v) < 2) return(helpText("patient has too few visits"))
    sliderInput("lm", "use history up to time", min(v), max(v),
                v[max(1, floor(length(v) / 3))], step = signif(diff(range(v)) / 20, 2))
  })

  pr <- reactive({
    req(rv$big, input$pid, input$lm)
    sid <- as.numeric(input$pid); L <- input$lm
    full <- subset(rv$dl, id == sid); fut <- subset(full, visit > L); nd <- nd_of(rv$dl, sid, L)
    PE <- predJMbayes(rv$big, ids = sid, newdata = nd, process = "event",
                      times = seq(L, rv$tmax, length.out = 30))
    ft <- sort(unique(c(seq(L, rv$tmax, length.out = 25), fut$visit))); ft <- ft[ft > L]
    ndL <- rbind(nd[, c("id", "visit", "x1", "y")],
                 data.frame(id = sid, visit = ft, x1 = nd$x1[1], y = NA_real_))
    ndL$time <- L; ndL$status <- 0
    PL <- predJMbayes(rv$big, ids = sid, newdata = ndL, process = "longitudinal")
    if (nrow(fut)) { r <- fut$y - patn(PL$p1[[1]], fut$visit); ty <- "held-out (forecast vs future)" }
    else { h <- subset(full, visit <= L); r <- h$y - patn(PL$p1[[1]], h$visit); ty <- "in-sample (fitted vs observed)" }
    list(PL = PL, PE = PE, full = full, fut = fut, L = L,
         rmse = sqrt(mean(r^2)), mae = mean(abs(r)), ty = ty, n = length(r))
  })

  output$pl <- renderPlot({
    if (is.null(rv$big)) { plot.new(); title("Generate data, then Fit model"); return() }
    p <- tryCatch(pr(), error = function(e) e)
    if (inherits(p, "error")) { plot.new(); text(.5, .5, paste(strwrap(conditionMessage(p), 55), collapse = "\n")); return() }
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)); on.exit(par(op))
    # left: biomarker history + forecast + dashed 95% CI
    fc <- as.data.frame(p$PL$p1[[1]]); o <- order(fc$visit)
    pc <- grep("^pred_", names(fc), value = TRUE)[1]
    lo <- grep("^low_", names(fc), value = TRUE)[1]; up <- grep("^upp_", names(fc), value = TRUE)[1]
    ho <- subset(p$full, visit <= p$L)
    plot(fc$visit[o], fc[[pc]][o], type = "n",
         ylim = range(c(fc[[lo]], fc[[up]], p$full$y), na.rm = TRUE),
         xlab = "time", ylab = "biomarker", main = "History + forecast")
    if (!is.na(lo)) {
      polygon(c(fc$visit[o], rev(fc$visit[o])), c(fc[[lo]][o], rev(fc[[up]][o])),
              col = "#2166ac22", border = NA)
      lines(fc$visit[o], fc[[lo]][o], col = "#2166ac", lty = 2)
      lines(fc$visit[o], fc[[up]][o], col = "#2166ac", lty = 2)
    }
    lines(fc$visit[o], fc[[pc]][o], col = "#2166ac", lwd = 2)
    points(ho$visit, ho$y, pch = 16)
    if (nrow(p$fut)) points(p$fut$visit, p$fut$y, pch = 17, col = "#b2182b")
    abline(v = p$L, lty = 2, col = "grey50")
    # right: conditional survival + dashed 95% CI
    ev <- as.data.frame(p$PE$p1[[1]]); Tv <- attr(p$PE$p1[[1]], "Time_var")
    if (is.null(Tv) || !Tv %in% names(ev)) Tv <- "time"
    o2 <- order(ev[[Tv]]); tt <- ev[[Tv]][o2]
    s <- 1 - ev$pred_CIF[o2]; slo <- 1 - ev$upp_CIF[o2]; sup <- 1 - ev$low_CIF[o2]
    plot(tt, s, type = "n", ylim = c(0, 1),
         xlab = "time", ylab = "survival (1 - CIF)", main = "Conditional survival")
    polygon(c(tt, rev(tt)), c(slo, rev(sup)), col = "#b2182b22", border = NA)
    lines(tt, slo, col = "#b2182b", lty = 2); lines(tt, sup, col = "#b2182b", lty = 2)
    lines(tt, s, col = "#b2182b", lwd = 2)
    abline(v = p$L, lty = 2, col = "grey50")
  })

  output$err <- renderText({
    if (is.null(rv$big)) return("Fit a model first.")
    p <- tryCatch(pr(), error = function(e) NULL); req(p)
    sprintf("error type   : %s\npoints scored: %d\nRMSE : %.3f\nMAE  : %.3f",
            p$ty, p$n, p$rmse, p$mae)
  })
}

shinyApp(ui, server)
