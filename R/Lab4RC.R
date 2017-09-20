linreg <- setRefClass("linreg",
                        fields = list(formula = "formula", data = "data.frame",Coefficients = "matrix", Fits = "matrix",
                                      Residual = "matrix", degreesFreedom = "numeric",
                                      rvariance = "matrix",
                                      var_betas="matrix",tB="matrix",DataName="character"),
                        
                        methods = list(
                          initialize = function(formula = formula, data = data){
                            mx <- model.matrix(formula, data=data)
                            y <- all.vars(formula)[1]
                            
                            ta <- t(mx)
                            
                            # Regressions coefficients
                            betas <- solve((ta %*% mx)) %*% ta %*% data[,y[1]==names(data)]
                            Coefficients <<- betas
                            
                            # The fitted values
                            y_hat <- mx %*% betas
                            Fits <<- y_hat
                            
                            # The residuals
                            e_hat <- data[,y[1]==names(data)] - y_hat
                            Residual <<- e_hat
                            
                            # The degrees of freedom
                            n <- nrow(mx)
                            p <- ncol(mx)
                            df <- n - p
                            degreesFreedom <<- df
                            
                            # The residual variance
                            rvariance <<- t(e_hat) %*% e_hat / df
                            
                            # The variance of the regression coefficients:
                            right_side <- 1/(t(mx)%*%mx)
                            var_betas <<- rvariance[1,1] * right_side
                            
                            # The t-values for each coefficient
                            tB <<- betas / sqrt(diag(var_betas))
                            formula <<- formula
                            DataName<<- deparse(substitute(data))
                          },
                          print = function(){
                            cat("Call:",sep="\n")
                            cat(paste("linreg(","formula = ",formula[2]," ",formula[1]," ",formula[3],", ","data = ",DataName,")",sep=""), sep="\n")
                            cat(sep="\n")
                            cat("Coefficients:")
                            cat(sep="\n")

                            namn<-rownames(Coefficients)
                            cat(" ")
                            cat(namn)
                            cat(" ")
                            cat(sep="\n")

                            cat(Coefficients)
                            
                            # return(
                            #   list(
                            #     Formula_call = formula,
                            #     Regression_Coefficient=Coefficients
                            #   )
                            # )
                            
                          },
                          plot = function(){
                            require(ggplot2)
                            
                            dataint <- data.frame(residual = Residual, fitos = Fits, std_residual = sqrt(abs(scale(Residual))))
                            ggplot(data = dataint, aes(x = fitos, y = residual) ) +
                              geom_point() + labs(x = "Fitted values", y = "Residuals") +
                              geom_smooth(method="loess", se = FALSE, color = "red") +
                              geom_hline(yintercept = 0) + theme_bw() + ggtitle("Residuals vs Fitted") +
                              theme(plot.title = element_text(hjust = 0.5))
                          },
                          resid = function(){
                            return(
                              Residual
                            )
                          },
                          pred = function(){
                            return(
                              Fits
                            )
                          },
                          coef = function(){
                            return(
                              Coefficients
                            )
                          }
                        ))

linreg_mod <- linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)
linreg_mod$print()
# linreg_mod$plot()
# linreg_mod$resid()
# linreg_mod$pred()
linreg_mod$coef()

