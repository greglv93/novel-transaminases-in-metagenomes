#setwd([insert wd here])
dops_incthresh = read.table('./dops_vs_incthresh_3.40.640.10.csv', sep=',', header=TRUE)

plot(dops_incthresh$dops, dops_incthresh$eval, log='y', ylim = c(1e-220, 1e-35),
     xlab = 'DOPS', ylab = 'log10(inclusion threshold)', cex=0.5)

log10_eval = log10(dops_incthresh$eval)
abline(lm(log10_eval~dops_incthresh$dops), col='red')
abline(v=70, lty=2)
cor.test(log10_eval, dops_incthresh$dops, method='pearson')
text(90, 1e-180, 'r=0.81', cex=0.8)
text(92.5, 1e-200, 'p=2.2e-16', cex=0.8)
legend('bottom', 'DOPS=70', lty=2, cex=0.8)

# length vs inc_thresh  for the 16 funfams with dops=0
  # (to check if length could be added into a multiple linear regression model)

len_incthresh = read.table('./lenvsincthresh_dops0.csv', sep=',', header = TRUE)
plot(len_incthresh$length, len_incthresh$eval, log='y', xlab = 'funfam length', 
     ylab = 'inclusion threshold')
cor.test(log10(len_incthresh$eval), len_incthresh$length, method='pearson')

######

ff_table = read.table('./v4.1_3.40.640.10_ffs_dops_eval_len.csv', sep=',', header = TRUE)
negativelog10_eval = -log10_eval
negativelog10eval_over_len = negativelog10_eval/(ff_table$len)
len_over_negativelog10eval = (ff_table$len)/negativelog10_eval
log10eval_over_len = log10_eval/(ff_table$len)
len_over_log10eval = ff(table$len)/log10_eval

######

plot(ff_table$dops, negativelog10eval_over_len,
     xlab = 'DOPS', ylab = '-log10(inc_thresh) / domain length', cex=0.5)
abline(v=70, lty=2)
abline(h=0.45, lty=2)
abline(lm(negativelog10eval_over_len~ff_table$dops), col='red')
cor.test(negativelog10eval_over_len, ff_table$dops, method='pearson')
text(30, 0.25, 'r=-0.94', cex=0.8)
text(32.2, 0.2, 'p=2.2e-16', cex=0.8)

######

plot(ff_table$dops, len_over_negativelog10eval,
     xlab = 'DOPS', ylab = 'domain length / -log10(inc_thresh)', cex=0.5)
abline(v=70, lty=2)
abline(h=2.2, lty=2)
abline(lm(len_over_negativelog10eval~ff_table$dops), col='red')
cor.test(ff_table$dops, len_over_negativelog10eval, method = 'pearson')

#######

plot(ff_table$dops, log10eval_over_len,
     xlab = 'DOPS', ylab = 'log10(inc_thresh) / domain length', cex=0.5)
abline(v=70, lty=2)
abline(h=-0.46, lty=2)
abline(lm(log10eval_over_len~ff_table$dops), col='red')
cor.test(log10eval_over_len, ff_table$dops, method='pearson')
text(30, -0.2, 'r=0.94', cex=0.8)
text(30, -0.25, 'p=2.2e-16', cex=0.8)
