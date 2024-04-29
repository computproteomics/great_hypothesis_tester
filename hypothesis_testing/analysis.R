
# ------------------------------------------------------- #
# Setup                                                   #
# ------------------------------------------------------- #

# get args from command
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop("Error: Expected exactly 2 args!")
}

message(paste("Starting Retention Time Analysis for", args[1], "..."))

# ------------------------------------------------------- #
# Analysis of retention time results                      #
# ------------------------------------------------------- #

# read and scale data
rtData <- read.csv(args[1])

prositRtData <- rtData$Prosit.Prediction
deeplcRtData <- rtData$DeepLC.Prediction

scaledPrositRtData <- (prositRtData - min(prositRtData)) / (max(prositRtData) - min(prositRtData))
scaledDeeplcRtData <- (deeplcRtData - min(deeplcRtData)) / (max(deeplcRtData) - min(deeplcRtData))

# plot histograms
color <- rgb(173, 216, 230, max=255, alpha=80, names='lt.blue')
ax <- pretty(0:1, n=12)

rt_histProsit <- hist(scaledPrositRtData, breaks=ax, plot=FALSE)
rt_histDeeplc <- hist(scaledDeeplcRtData, breaks=ax, plot=FALSE)

plot(rt_histProsit, col=color, xlab='Retention Time', ylab='Number of Listings', main='Prosit RT predictions')
plot(rt_histDeeplc, col=color, xlab='Retention Time', ylab='Number of Listings', main='DeepLC RT predictions')

# plot boxplots
rt_boxplot_data <- list(
    Prosit = scaledPrositRtData,
    DeepLC = scaledDeeplcRtData
)

rt_boxplot = boxplot(rt_boxplot_data, col='lightblue', main='Retention Time per Model', xlab='Model', ylab='Retention Time')

# ------------------------------------------------------- #
# Analysis of fragmentation  results                      #
# ------------------------------------------------------- #

# read fragmentation pattern data
fragData <- read.csv('C:/Users/user/Desktop/Studium/Berufspraktikum/great_hypothesis_tester/results/100k_10_length_random_sequences_results/100k_10_length_random_sequences_fragmentation_results.csv')

# filter m/z and intensity values
yMzProsit <- fragData[, grep('y\\d+_mz_prosit', names(fragData))]
bMzProsit <- fragData[, grep('b\\d+_mz_prosit', names(fragData))]

yMzDeeplc <- fragData[, grep('y\\d+_mz_deeplc', names(fragData))]
bMzDeeplc <- fragData[, grep('b\\d+_mz_deeplc', names(fragData))]

yIntensProsit <- fragData[, grep('y\\d+_intens_prosit', names(fragData))]
bIntensProsit <- fragData[, grep('b\\d+_intens_prosit', names(fragData))]

yIntensDeeplc <- fragData[, grep('y\\d+_intens_deeplc', names(fragData))]
bIntensDeeplc <- fragData[, grep('b\\d+_intens_deeplc', names(fragData))]

# setup for histograms
figures_per_row <- 3
amount_rows <- ceiling(length(yMzProsit) / figures_per_row)

# setup for boxplots
y_names <- paste0('y', seq_along(names(yMzProsit)))
b_names <- paste0('b', seq_along(names(yMzProsit)))

names(yMzProsit) <- y_names
names(yMzDeeplc) <- y_names
names(bMzProsit) <- b_names
names(bMzDeeplc) <- b_names

# print histograms for m/z values 
par(mfrow=c(amount_rows, figures_per_row))

# y-ions prosit
for (i in colnames(yMzProsit)) {
  hist(yMzProsit[[i]], xlab='M/z value', ylab='Frequency',
       main=paste(strsplit(i, split='_')[[1]][1], 'm/z predictions Prosit'))
}

# b-ions prosit
for (i in colnames(bMzProsit)) {
  hist(bMzProsit[[i]], xlab='m/z value', ylab='Frequency', 
        main=paste(strsplit(i, split='_')[[1]][1], 'm/z predictions Prosit'))
}

# y-ions deeplc
for (i in colnames(yMzDeeplc)) {
  hist(yMzDeeplc[[i]], xlab='m/z value', ylab='Frequency', 
        main=paste(strsplit(i, split='_')[[1]][1], 'm/z predictions DeepLC'))
}

# b-ions deeplc
for (i in colnames(bMzDeeplc)) {
  hist(bMzDeeplc[[i]], xlab='m/z value', ylab='Frequency', 
        main=paste(strsplit(i, split='_')[[1]][1], 'm/z predictions DeepLC'))
}

# print boxplots for m/z values
par(mfrow=c(1,1))

# y-ions prosit
boxplot(yMzProsit, col=rainbow(ncol(yMzProsit)), main='m/z Prosit predictions for y-ions')

# b-ions prosit
boxplot(bMzProsit, col=rainbow(ncol(bMzProsit)), main='m/z Prosit predictions for b-ions')

# y-ions deeplc
boxplot(yMzDeeplc, col=rainbow(ncol(yMzDeeplc)), main='m/z DeepLC predictions for y-ions')

# b-ions deeplc
boxplot(bMzDeeplc, col=rainbow(ncol(bMzDeeplc)), main='m/z DeepLC predictions for b-ions')

# print histograms for intensity values
par(mfrow=c(amount_rows, figures_per_row))

# y-ions prosit
for (i in colnames(yIntensProsit)) {
  hist(yIntensProsit[[i]], xlab='Intensity', ylab='Frequency',
       main=paste(strsplit(i, split='_')[[1]][1], 'intensity predictions Prosit'))
}

# b-ions prosit
for (i in colnames(bIntensProsit)) {
  hist(bIntensProsit[[i]], xlab='Intensity', ylab='Frequency', 
        main=paste(strsplit(i, split='_')[[1]][1], 'intensity predictions Prosit'))
}

# y-ions deeplc
for (i in colnames(yIntensDeeplc)) {
  hist(yIntensDeeplc[[i]], xlab='Intensity', ylab='Frequency', 
        main=paste(strsplit(i, split='_')[[1]][1], 'intensity predictions DeepLC'))
}

# b-ions deeplc
for (i in colnames(bIntensDeeplc)) {
  hist(bIntensDeeplc[[i]], xlab='Intensity', ylab='Frequency', 
        main=paste(strsplit(i, split='_')[[1]][1], 'intensity predictions DeepLC'))
}

# print boxplots for intensity values
par(mfrow=c(1,1))

# y-ions prosit
boxplot(yIntensProsit, col=rainbow(ncol(yIntensProsit)), main='Intensity Prosit predictions for y-ions')

# b-ions prosit
boxplot(bIntensProsit, col=rainbow(ncol(bIntensProsit)), main='Intensity Prosit predictions for b-ions')

# y-ions deeplc
boxplot(yIntensDeeplc, col=rainbow(ncol(yIntensDeeplc)), main='Intensity DeepLC predictions for y-ions')

# b-ions deeplc
boxplot(bIntensDeeplc, col=rainbow(ncol(bIntensDeeplc)), main='Intensity DeepLC predictions for b-ions')
