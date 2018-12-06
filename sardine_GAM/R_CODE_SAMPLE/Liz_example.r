library(ROCR)
library(mgcv)


load('revisedspawninghab.rdata')


# logistic gam, limit splines to unimodel curves with k=3 knots,
# 'select' and 'gamma' args are for penalized smoothing -- predictors
# can be shrunken to zero and removed from the model
mod <- gam(sardpres ~ s(chlamn50, bs='cs', k=3) + 
           s(chlamaxdepth, bs='cs', k=3) + 
           s(tempmn50, bs='cs', k=3) + 
           s(salinmn50, bs='cs', k=3) + 
           s(dayofyear, bs='cs', k=3) + 
           s(slope, bs='cs', k=3) + 
           sardprevdens, scale=0, 
           data=dat, family='binomial', gamma=1.4, 
           select=TRUE) 
summary(mod)
par(mfrow=c(2, 4))
plot(mod, se=2, all.terms=TRUE, trans=plogis)

# optionally refit a final model with predictors shrunken to linear changed
# to linear predictors instead of splines, and predictors shrunken to
# zero removed
mod <- gam(sardpres ~ s(chlamn50, bs='cs', k=3) + 
           tempmn50 + salinmn50 + sardprevdens, scale=0, 
           data=dat, family='binomial', gamma=1.4, 
           select=TRUE) 
summary(mod)
par(mfrow=c(1, 4))
plot(mod, se=2, all.terms=TRUE, trans=plogis)


# get modeled (predicted) values for the original data as proportions
plogis(predict(mod, dat))


# I can no longer remember why I did not allow interactions in that
# paper, but you can also do something like this:
mod2 <- gam(sardpres ~ te(tempmn50, chlamn50, bs=c('cs', 'cs')) +
            salinmn50 + sardprevdens, scale=0, 
            data=dat, family='binomial', gamma=1.4, 
            select=TRUE)
summary(mod2)
plot(mod2, residuals=T)

# looks like the original model was better in this case though:
AIC(mod, mod2)

