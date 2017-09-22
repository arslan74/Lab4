#' @title Linear Regression
#' @description You can have Reference Class containing some calculations by giving formula and data.
#' @field formula Formula
#' @field data A data frame
#' @examples
#' data(iris)
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$print()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$pred()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$summary()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$resid()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$coef()
#' linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)$plot()
#' @export linreg
linreg <- setRefClass("linreg",
                        fields = list(formula = "formula", data = "data.frame",Coefficients = "numeric", Fits = "numeric",
                                      Residuals = "numeric", df = "numeric",
                                      rvariance = "matrix",
                                      Var_residuals="numeric",tBetas="numeric",DataName="character",Std_betas="numeric",Pvalues="numeric"),
                        
                        methods = list(
                          initialize = function(formula = formula, data = data){
                            
                            formula<<-formula
                            data<<-data

                            x<-model.matrix(formula,data)
                            y<-all.vars(formula)[1]
                            y<-as.matrix(data[,names(data)==y])
                            
                            b_hat<-solve(t(x)%*%x)%*%t(x)%*%y
                            y_fits<-x%*%b_hat
                            e<-y-y_fits
                            df1<-length(y)-ncol(x)
                            var_e<-(t(e)%*%e)/df1
                            
                            var_b_hat<-as.numeric(var_e)*diag(solve((t(x)%*%x)))
                            std_b_hat<-sqrt(var_b_hat)
                            t_b_hat<-as.numeric(b_hat)/std_b_hat
                            p_b_hat<-(1-(pt(abs(t_b_hat),df = df1)))*2
                            
                            b_hat_numeric<-as.numeric(b_hat)
                            names(b_hat_numeric)<-rownames(b_hat)
                            
                            input_var<-as.character(match.call(expand.dots = FALSE))
                            Coefficients<<-b_hat_numeric
                            Fits<<-as.numeric(y_fits)
                            Residuals<<-as.numeric(e)
                            df<<-df1
                            Var_residuals<<-as.numeric(var_e)
                            Std_betas<<-std_b_hat
                            tBetas<<-t_b_hat
                            Pvalues<<-p_b_hat
                            DataName<<- deparse(substitute(data))
                            
                          },
                          #' @title Print regression coefficients
                          #' @name  print
                          #' @description This function prints regression coefficients by using given formula and data in initialization.
                          print = function(){
                            cat("Call:",sep="\n")
                            cat(paste("linreg(","formula = ",formula[2]," ",formula[1]," ",formula[3],", ","data = ",DataName,")",sep=""), sep="\n")
                            cat(sep="\n")
                            cat("Coefficients:")
                            cat(sep="\n")

                            beta<-Coefficients
                            namn<-names(beta)
                            names(beta)<-NULL
                            beta<-round(beta,4)
                            
                            for(i in 2:length(beta)){
                              beta[i]<-format(beta[i], width=max(nchar(beta[i]),nchar(namn[i])),justify = "right")
                            }

                            beta[1]<-format(beta[1], width=max(nchar(beta[1]),nchar(namn[1]),nchar("Coefficients")),justify = "right")
                            namn[1]<-format(namn[1], width=max(nchar(beta[1]),nchar(namn[1]),nchar("Coefficients")),justify = "right")

                            beta[1]<-paste(beta[1],"  ",sep="")
                            namn[1]<-paste(namn[1],"  ",sep="")

                            beta[2]<-paste(beta[2]," ",sep="")
                            namn[2]<-paste(namn[2]," ",sep="")

                            cat(" ")
                            cat(namn)
                            cat(" ")
                            cat(sep="\n")
                            cat(beta)
                            
                          },
                          #' @title Plot graphs
                          #' @name  plot
                          #' @description This function plots two graphs, such as Fitted values vs Residuals and Scale Location by using given formula and data in initialization.
                          plot = function(){
                            require(ggplot2)
                            dataint <- data.frame(residual = Residuals, fitos = Fits)
                            a <- ggplot(data = dataint, aes(x = fitos, y = residual) ) +
                              geom_point() + labs(x = "Fitted values", y = "Residuals") +
                              geom_smooth(method="loess", se = FALSE, color = "red") +
                              geom_hline(yintercept = 0) + theme_bw() + ggtitle("Residuals vs Fitted") +
                              theme(plot.title = element_text(hjust = 0.5))


                            dataint2 <- data.frame(residual = sqrt(abs(Residuals)), fitos = Fits)
                            b <- ggplot(data = dataint2, aes(x = fitos, y = residual) ) +
                              geom_point() + labs(x = "Fitted values", y = expression(sqrt(abs("Standardized residuals")))) +
                              geom_smooth(method="loess", se = FALSE, color = "red") +
                              geom_hline(yintercept = 0) + theme_bw() + ggtitle("Scale Location") +
                              theme(plot.title = element_text(hjust = 0.5))
                            
                            return(list(ResidualsVsFitted = a, ScaleLocation = b))
                            
                          },
                          #' @title Residuals
                          #' @name  resid
                          #' @return Residuals
                          #' @description This function returns residuals value.
                          resid = function(){
                            return(
                              Residuals
                            )
                          },
                          #' @title Fitted Values
                          #' @name  pred
                          #' @return Fits
                          #' @description This function returns fitted value.
                          pred = function(){
                            return(
                              Fits
                            )
                          },
                          #' @title Regressions Coefficients
                          #' @name  coef
                          #' @return Coefficients
                          #' @description This function returns regression coefficients
                          coef = function(){
                            return(
                              Coefficients
                            )
                          },
                          #' @title Summary
                          #' @name  summary
                          #' @description This function prints the coefficients with their standard error, t-value and p-value.
                          summary = function(){
                            # (Intercept) -2.55 0.55 -4.44 ****
                            # Sepal.Width -1.32 0.15 -10.95  ****
                            # Sepal.Length 1.72 0.05 27.55  ****
                            # Residual standard error: 0.63 on 147 degrees of freedom
                            
                            beta<-Coefficients
                            namn<-names(beta)
                            names(beta)<-NULL
                            beta<-round(beta,4)
                            
                            for(i in 1:length(beta)){
                              beta[i]<-format(beta[i], width=max(nchar(beta[i]),nchar(namn[i])),justify = c("right"))
                              namn[i]<-format(namn[i], width=max(nchar(beta[i]),nchar(namn[i])),justify = c("right"))
                            }
                            
                            Variable<-as.character(names(Coefficients))
                            Estimate<-round(Coefficients,3)
                            Std_Error<-round(Std_betas,3)
                            t_value<-round(tBetas,3)
                            Pr_t<-round(Pvalues,5)
                            
                            svar<-data.frame(Variable,Estimate,Std_Error,t_value,Pr_t)
                            row.names(svar)<-NULL
                            svar$Variable<-as.character(svar$Variable)

                            cat("Call:",sep="\n")
                            cat(paste("linreg(","formula = ",formula[2]," ",formula[1]," ",formula[3],", ","data = ",DataName,")",sep=""), sep="\n")
                            cat(sep="\n")
                            cat("Coefficients:",sep="\n")
                            cat()
                            for(i in 1:nrow(svar)){
                              cat(paste(svar[i,],collapse = " "),sep="",collapse=" ***\n")
                            }
                            cat("",sep="\n")
                            cat(paste("Residual standard error: ",round(sqrt(Var_residuals),5) ," on " ,df, " degrees of freedom",sep=""))
                          } 
                        ))

linreg_mod <- linreg$new(Petal.Length~Species, data=iris)

# linreg_mod$print()
# linreg_mod$plot()
# linreg_mod$resid()
# linreg_mod$pred()
# linreg_mod$coef()
# linreg_mod$summary()