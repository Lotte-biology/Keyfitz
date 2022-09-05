## RAW CONTRAST - COMADRE -----------------------

plot(comadreEntropy$NewEntropy,
     comadreEntropy$OriginalEntropy,
     pch = 16, col = alpha("white", 0.3), cex = 1,
     ylim = c(0, 2.1), xlim = c(0, 2.1),
     ylab = "Original Entropy", xlab = "New Entropy")


rect(-1,1,1,3,
     col = alpha("blue", 0.05))
rect(1,-1,3,1,
     col = alpha("blue", 0.05))
rect(-1,-1,1,1,
     col = alpha("grey", 0.05))
rect(1,1,3,3,
     col = alpha("grey", 0.05))


abline(h = 1, col = "white")
abline(v = 1, col = "white")

criticalZone <- which(comadreEntropy$NewEntropy>1 & comadreEntropy$Original<1)

dev <- comadreEntropy$OriginalEntropy-comadreEntropy$NewEntropy
for(i in 1:length(comadreEntropy$NewEntropy)){
  lines(c(comadreEntropy$NewEntropy[[i]],
          comadreEntropy$NewEntropy[[i]]),
        c(comadreEntropy$NewEntropy[[i]], 
          comadreEntropy$NewEntropy[[i]]+dev[[i]]),
        col = alpha("black", 0.1))
}

points(comadreEntropy$NewEntropy[-criticalZone],
     comadreEntropy$OriginalEntropy[-criticalZone],
     pch = 16, col = alpha("blue", 0.2), cex = 1,
     ylim = c(0, 2.3), xlim = c(0, 2.3))

points(comadreEntropy$NewEntropy[-criticalZone],
       comadreEntropy$OriginalEntropy[-criticalZone],
       col = alpha("blue", 0.5), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(seq(-1,3),
       seq(-1,3),
       lty = 2, type = "l", col = "blue")

points(comadreEntropy$NewEntropy[criticalZone],
       comadreEntropy$OriginalEntropy[criticalZone],
       col = alpha("darkblue", 0.6),
       xlim = c(0,2), ylim = c(0,2), pch = 16, cex = 1)

points(comadreEntropy$NewEntropy[criticalZone],
       comadreEntropy$OriginalEntropy[criticalZone],
       col = alpha("darkblue", 0.5),
       xlim = c(0,2), ylim = c(0,2), cex = 1)


points(stageEntropyComadre$NewEntropy,
       stageEntropyComadre$OriginalEntropy,
       col = alpha("brown1", 0.8),
       xlim = c(0,2), ylim = c(0,2), pch=16, cex = 1)

points(stageEntropyComadre$NewEntropy,
       stageEntropyComadre$OriginalEntropy,
       col = alpha("brown1"),
       xlim = c(0,2), ylim = c(0,2), cex = 1)



Fit <- lm(OriginalEntropy~exp(NewEntropy), data = comadreEntropy)
summary(Fit)
# abline(Fit,lwd=2,col=alpha("red", 0.5))

newx <- seq(-1, 3, by = 0.01)
conf_interval <- predict(Fit,
                         newdata=data.frame(NewEntropy = newx),
                         interval="confidence",
                         level = 0.95)

points(newx, conf_interval[,1], type = "l", col = alpha("red", 0.4))

polygon(c(newx, rev(newx)),
        c(conf_interval[,2],
          rev(conf_interval[,3])),
        col = alpha("red", 0.2),
        border = alpha("red", 0.1))
##-----------

plot(comadreEntropy$NewEntropy[-8], dev[-8],
     cex = 0.5,
     col = "white",
     xlab = "New Entropy",
     ylab = "Distance to Original Entropy",
     xlim = c(0, 2.1))

for(i in 1:length(comadreEntropy$NewEntropy)){
  lines(c(comadreEntropy$NewEntropy[[i]],
          comadreEntropy$NewEntropy[[i]]),
        c(dev[[i]], 
          0),
        col = alpha("blue", 0.1),
        lwd = 2)
}

abline(h = 0, col = "blue", lty = 2)


##---------

gg_comadreEntropy <- as.data.frame(comadreEntropy)
gg_comadreEntropy <- gg_comadreEntropy[-8,]


base1 <- ggplot(gg_comadreEntropy, aes(x=NewEntropy, y=OriginalEntropy) ) +
  xlim(0, 2.1)+
  ylim(0, 2.1)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  theme(legend.position='none')+
  theme_minimal()+
  scale_fill_gradient(
    low = alpha("white", 0), 
    high = "blue"
  )

base2 <- base1 + geom_point(alpha = 0.1 )

base3 <- base2 +
  stat_smooth(method = "lm",
              formula = y ~ exp(x),
              geom = "smooth",
              color = "brown1",
              se = FALSE)

base3 +
  geom_hline(yintercept=1, linetype="dashed",
             color = "black", size=0.5) +
  geom_vline(xintercept=1, linetype="dashed",
                   color = "black", size=0.5)+
  geom_abline(intercept = 0, slope = 1, color="blue", 
                linetype="dashed", size=0.5)




## RAW CONTRAST - COMPADRE ------------------------------------------------------------

plot(compadreEntropy$NewEntropy,
     compadreEntropy$OriginalEntropy,
     pch = 16, col = alpha("white", 0.3), cex = 1,
     ylim = c(0, 2.1), xlim = c(0, 2.1),
     ylab = "Original Entropy", xlab = "New Entropy")


rect(-1,1,1,3,
     col = alpha("blue", 0.05))
rect(1,-1,3,1,
     col = alpha("blue", 0.05))
rect(-1,-1,1,1,
     col = alpha("grey", 0.05))
rect(1,1,3,3,
     col = alpha("grey", 0.05))

abline(h = 1, col = "white")
abline(v = 1, col = "white")

criticalZone <- which(compadreEntropy$NewEntropy>1 & compadreEntropy$Original<1)

dev <- compadreEntropy$OriginalEntropy-compadreEntropy$NewEntropy
for(i in 1:length(compadreEntropy$NewEntropy)){
  lines(c(compadreEntropy$NewEntropy[[i]],
          compadreEntropy$NewEntropy[[i]]),
        c(compadreEntropy$NewEntropy[[i]], 
          compadreEntropy$NewEntropy[[i]]+dev[[i]]),
        col = alpha("black", 0.1))
}

points(compadreEntropy$NewEntropy[-criticalZone],
       compadreEntropy$OriginalEntropy[-criticalZone],
       pch = 16, col = alpha("blue", 0.2), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(compadreEntropy$NewEntropy[-criticalZone],
       compadreEntropy$OriginalEntropy[-criticalZone],
       col = alpha("blue", 0.5), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(seq(-1,3),
       seq(-1,3),
       lty = 2, type = "l", col = "blue")

points(compadreEntropy$NewEntropy[criticalZone],
       compadreEntropy$OriginalEntropy[criticalZone],
       col = alpha("darkblue", 0.6),
       xlim = c(0,2), ylim = c(0,2), pch = 16, cex = 1)

points(compadreEntropy$NewEntropy[criticalZone],
       compadreEntropy$OriginalEntropy[criticalZone],
       col = alpha("darkblue", 0.5),
       xlim = c(0,2), ylim = c(0,2), cex = 1)


points(stageEntropyCompadre$NewEntropy,
       stageEntropyCompadre$OriginalEntropy,
       col = alpha("brown1", 0.8),
       xlim = c(0,2), ylim = c(0,2), pch=16, cex = 1)

points(stageEntropyCompadre$NewEntropy,
       stageEntropyCompadre$OriginalEntropy,
       col = alpha("brown1"),
       xlim = c(0,2), ylim = c(0,2), cex = 1)


backgroundEntropyCompadre <- calculateEntropy(compadreMaxDim)

points(backgroundEntropyCompadre$NewEntropy,
       backgroundEntropyCompadre$OriginalEntropy,
       col = alpha("gray", 0.7),
       xlim = c(0,2), ylim = c(0,2), pch=16, cex = 1)

points(backgroundEntropyCompadre$NewEntropy,
       backgroundEntropyCompadre$OriginalEntropy,
       col = alpha("gray", 0.7),
       xlim = c(0,2), ylim = c(0,2), cex = 1)


Fit <- lm(OriginalEntropy~NewEntropy, data = stageEntropyCompadre)
summary(Fit)
abline(Fit,lwd=2,col=alpha("red", 0.5))

newx <- seq(-1, 3, by = 0.01)
conf_interval <- predict(Fit,
                         newdata=data.frame(NewEntropy = newx),
                         interval="confidence",
                         level = 0.95)

polygon(c(newx, rev(newx)),
        c(conf_interval[,2],
          rev(conf_interval[,3])),
        col = alpha("red", 0.2),
        border = alpha("red", 0.1))



Fit <- lm(OriginalEntropy~NewEntropy, data = backgroundEntropyCompadre)
summary(Fit)
abline(Fit,lwd=2,col=alpha("gray", 0.5))

newx <- seq(-1, 3, by = 0.01)
conf_interval <- predict(Fit,
                         newdata=data.frame(NewEntropy = newx),
                         interval="confidence",
                         level = 0.95)

polygon(c(newx, rev(newx)),
        c(conf_interval[,2],
          rev(conf_interval[,3])),
        col = alpha("gray", 0.2),
        border = alpha("gray", 0.1))



# 



## MEASURE CONTRAST - COMADRE ------------------------------------------------------------


criticalZone <- which(comadreEntropy$NewEntropy>1 & comadreEntropy$Original<1)

plot(log(comadreEntropy$lifeExpOut),
     comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy,
     pch = 16, col = alpha("white", 0.3), cex = 1,
     ylab = "Different in Entropy Measures", xlab = "Log Life Expectancy",
     ylim = c(-0.75, 1), xlim = c(0, 5))

points(log(comadreEntropy$lifeExpOut)[-criticalZone],
       (comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy)[-criticalZone],
       pch = 16, col = alpha("blue", 0.2), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(log(comadreEntropy$lifeExpOut)[-criticalZone],
       (comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy)[-criticalZone],
       col = alpha("blue", 0.2), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(log(comadreEntropy$lifeExpOut)[criticalZone],
       (comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy)[criticalZone],
       pch = 16, col = alpha("darkblue", 0.6), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(log(comadreEntropy$lifeExpOut)[criticalZone],
       (comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy)[criticalZone],
       col = alpha("darkblue", 0.9), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

abline(h = 0, col = "blue", lty = 2)


library(npreg)

# fit using ss

holdX <- log(comadreEntropy$lifeExpOut)[criticalZone]
holdY <- (comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy)[criticalZone]
mod.ss <- ss(holdX, holdY, nknots = 5)
lines(holdX, mod.ss$y, lty = 2, col = 2, lwd = 2)
lines(mod.ss$x, mod.ss$y, col = "red", lwd = 2)


# Non Logged -----------------------------------------------------------------------------

criticalZone <- which(comadreEntropy$NewEntropy>1 & comadreEntropy$Original<1)

plot(comadreEntropy$lifeExpOut,
     comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy,
     pch = 16, col = alpha("white", 0.3), cex = 1,
     ylab = "Different in Entropy Measures", xlab = "Log Life Expectancy",
     ylim = c(-0.75, 1), xlim = c(0, 20))

points(comadreEntropy$lifeExpOut[-criticalZone],
       (comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy)[-criticalZone],
       pch = 16, col = alpha("blue", 0.2), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(log(comadreEntropy$lifeExpOut)[-criticalZone],
       (comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy)[-criticalZone],
       col = alpha("blue", 0.2), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(log(comadreEntropy$lifeExpOut)[criticalZone],
       (comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy)[criticalZone],
       pch = 16, col = alpha("darkblue", 0.6), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(log(comadreEntropy$lifeExpOut)[criticalZone],
       (comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy)[criticalZone],
       col = alpha("darkblue", 0.9), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

abline(h = 0, col = "blue", lty = 2)







## MEASURE CONTRAST - COMPADRE ------------------------------------------------------------


criticalZone <- which(compadreEntropy$NewEntropy>1 & compadreEntropy$Original<1)

plot(log(compadreEntropy$lifeExpOut),
     compadreEntropy$NewEntropy-compadreEntropy$OriginalEntropy,
     pch = 16, col = alpha("white", 0.3), cex = 1,
     ylab = "Different in Entropy Measures", xlab = "Log Life Expectancy",
     ylim = c(-1.5, 1), xlim = c(0, 5))

points(log(compadreEntropy$lifeExpOut),
       (compadreEntropy$NewEntropy-compadreEntropy$OriginalEntropy),
       pch = 16, col = alpha("blue", 0.2), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

points(log(compadreEntropy$lifeExpOut),
       (compadreEntropy$NewEntropy-compadreEntropy$OriginalEntropy),
       col = alpha("blue", 0.2), cex = 1,
       ylim = c(0, 2.3), xlim = c(0, 2.3))

backgroundEntropyCompadre <- calculateEntropy(compadreMaxDim)

points(backgroundEntropyCompadre$lifeExpOut,
       backgroundEntropyCompadre$NewEntropy-backgroundEntropyCompadre$OriginalEntropy,
       col = alpha("gray", 0.7),
       xlim = c(0,2), ylim = c(0,2), pch=16, cex = 1)

points(backgroundEntropyCompadre$lifeExpOut,
       backgroundEntropyCompadre$NewEntropy-backgroundEntropyCompadre$OriginalEntropy,
       col = alpha("gray", 0.7),
       xlim = c(0,2), ylim = c(0,2),cex = 1)

abline(h = 0, col = "blue", lty = 2)





# # COMADRE Plot Difference Entropy Measures -----------
# 
# plot(comadreEntropy$lifeExpOut,
#      comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy,
#      pch=16, col = alpha("gray", 0.5), cex = 1,
#      ylim = c(-1, 1), xlim = c(1, 120))
# 
# points(comadreEntropy$lifeExpOut,
#        comadreEntropy$NewEntropy-comadreEntropy$OriginalEntropy,
#        col = alpha("gray", 0.5), cex = 1)
# 
# abline(h = 0, lty = 2, type = "l", col = "red")
# 
# 
# # COMPADRE Plot Comparative Entropy Measures -----------
# 
# plot(compadreEntropy$NewEntropy,
#      compadreEntropy$OriginalEntropy,
#      pch = 16, col = alpha("gray", 0.3), cex = 1,
#      ylim = c(0, 2.3), xlim = c(0, 2.3))
# 
# points(compadreEntropy$NewEntropy,
#        compadreEntropy$OriginalEntropy,
#        col = alpha("gray", 0.5), cex = 1,
#        ylim = c(0, 2.3), xlim = c(0, 2.3))
# 
# points(seq(-1,3),
#        seq(-1,3),
#        lty = 2, type = "l", col = "red")
# 
# # COMPADRE Plot Difference Entropy Measures -----------
# 
# plot(compadreEntropy$lifeExpOut,
#      evalCheck$compadreEntropy-compadreEntropy$OriginalEntropy,
#      pch=16, col = alpha("gray", 0.5), cex = 1,
#      ylim = c(-1, 1), xlim = c(1, 120))
# 
# points(compadreEntropy$lifeExpOut,
#        compadreEntropy$NewEntropy-compadreEntropy$OriginalEntropy,
#        col = alpha("gray", 0.5), cex = 1)
# 
# abline(h = 0, lty = 2, type = "l", col = "red")







