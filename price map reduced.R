library(GWmodel)
library(spdep)
library(sp)
library(tmap)
library(rgdal)
library(RColorBrewer)



## PRICE MAP
setwd("~/DTU/3. semestr/GIS/project/price map reduced")
# Load the output area shapefiles
pricemap<- readOGR(".", "price map reduced")

pricemap$CENA<-as.double(pricemap$CENA) 
pricemap$AREA_SQ_M<-as.double(pricemap$AREA_SQ_M) 

DeVar <- "CENA"
InDeVars <- c("AREA_SQ_M","HUBDISTBUS","HUBDISTMET","HUBDISTTRA","HUBDISTTR2")
model.sel <-model.selection.gwr(DeVar,InDeVars,data =pricemap, kernel="bisquare", adaptive=TRUE,bw=bw.gwr.1)

sorted.models <- model.sort.gwr(model.sel, numVars = length(InDeVars),ruler.vector = model.sel[[2]][,2])
model.list <- sorted.models[[1]]
model.view.gwr(DeVar, InDeVars, model.list = model.list)
plot(sorted.models[[2]][,2], col = "black", pch = 20, lty = 5,main = "Alternative view of GWR model selection procedure",
     ylab = "AICc", xlab = "Model number", type = "b")

#Basic GWR working
#Bisquare kernel
bw.gwr.1 <-bw.gwr(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                    HUBDISTTR2, data = pricemap,
                  approach = "AICc", kernel = "bisquare", adaptive = TRUE)

bw.gwr.cv <-bw.gwr(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                    HUBDISTTR2, data = pricemap,
                  approach = "CV", kernel = "bisquare", adaptive = TRUE)
bw.gwr.cv

gwr.res <- gwr.basic(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                     HUBDISTTR2, data = pricemap,
                     bw = bw.gwr.1,kernel = "bisquare", adaptive = TRUE, F123.test = FALSE)


print(gwr.res)
summary(gwr.res$SDF$Local_R2)


names(gwr.res$SDF)
mypalette.6 <- brewer.pal(6, "Spectral")
spplot(gwr.res$SDF, "AREA_SQ_M", key.space = "right",
         
            col.regions = mypalette.6, at = c(-8, -6, -4, -2, 0, 2, 4),
       
            main = "Basic GW regression coefficient estimates for AREA_SQ_M",
          
            sp.layout = map.layout)




print(gwr.res)
print(rgwr.res)


#Inspecting different kernels
#Gaussian
bw.gwr.gauss <-bw.gwr(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                        HUBDISTTR2, data = pricemap,
                      approach = "AICc", kernel = "gaussian", adaptive = TRUE)

gwr.gauss <- gwr.basic(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                         HUBDISTTR2, data = pricemap,
                       bw = bw.gwr.1,kernel = "gaussian", adaptive = TRUE, F123.test = FALSE)

summary(gwr.gauss$SDF$Local_R2)


#Exponential
bw.gwr.exp <-bw.gwr(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                      HUBDISTTR2, data = pricemap,
                    approach = "AICc", kernel = "exponential", adaptive = TRUE)

gwr.exp <- gwr.basic(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                       HUBDISTTR2, data = pricemap,
                     bw = bw.gwr.1,kernel = "exponential", adaptive = TRUE, F123.test = FALSE)

summary(gwr.exp$SDF$Local_R2)

#Tri-cube

bw.gwr.tricube <-bw.gwr(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                          HUBDISTTR2, data = pricemap,
                        approach = "AICc", kernel = "tricube", adaptive = TRUE)

gwr.tricube <- gwr.basic(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                           HUBDISTTR2, data = pricemap,
                         bw = bw.gwr.1,kernel = "tricube", adaptive = TRUE, F123.test = FALSE)

summary(gwr.tricube$SDF$Local_R2)


#bisquare has the lowest R2 median





#Inspecting colinearity
gwr.collin <- gwr.collin.diagno(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                                  HUBDISTTR2, data = pricemap, bw.gwr.1, kernel="bisquare",
                  adaptive=FALSE, p=2, theta=0, longlat=F)

summary(gwr.collin$local_CN)


#Robust GWR working
rgwr.res <- gwr.robust(CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                         HUBDISTTR2, data = pricemap,
                       bw = bw.gwr.1,kernel = "bisquare", adaptive = TRUE, F123.test = TRUE)


#Local compensated Ridge
library("car")
lm.global <- lm (CENA ~ AREA_SQ_M + HUBDISTBUS + HUBDISTMET + HUBDISTTRA +
                   HUBDISTTR2, data = pricemap)

vif(lm.global)
lcrm2.bw <-bw.gwr.lcr(CENA ~ AREA_SQ_M + HUBDISTBUS+ HUBDISTMET + HUBDISTTRA +
                          HUBDISTTR2 , data = pricemap,kernel = "bisquare", adaptive = TRUE)

lcrm2 <- gwr.lcr(CENA ~ AREA_SQ_M + HUBDISTBUS+ HUBDISTMET + HUBDISTTRA +
                     HUBDISTTR2 , data = pricemap, bw=lcrm2.bw, kernel="bisquare", adaptive=TRUE)


summary(gwr.lcr$SDF$Local_CN)
bw.gwr.cv


library("RColorBrewer")
map.na = list("SpatialPolygonsRescale", layout.north.arrow(),offset = c(329000, 261500), scale = 4000, col = 1)

map.scale.1 = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(326500, 217000), scale = 5000, col = 1,
                   fill = c("transparent", "blue"))
map.scale.2= list("sp.text", c(326500, 217900), "0", cex = 0.9, col = 1)
map.scale.3= list("sp.text", c(331500, 217900), "5km", cex = 0.9, col = 1)
map.layout <- list(map.na, map.scale.1, map.scale.2, map.scale.3)
mypalette.1 <- brewer.pal(8, "Reds")
mypalette.2 <- brewer.pal(5, "Blues")
mypalette.3 <- brewer.pal(6, "Greens")
mypalette.7 <- brewer.pal(8, "Reds")
spplot(gwr.lcr$SDF, "Local_CN", key.space = "right", col.regions = mypalette.7, at = seq(30, 110, length = 9),
            main = "Local condition numbers from basic GW regression",sp.layout = map.layout)


lcrm3.bw <- bw.gwr.lcr(CENA ~ AREA_SQ_M + HUBDISTBUS+ HUBDISTMET + HUBDISTTRA +
                   HUBDISTTR2 , data = pricemap, kernel = "bisquare", adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30)
lcrm3.bw

lcrm3 <- gwr.lcr(CENA ~ AREA_SQ_M + HUBDISTBUS+ HUBDISTMET + HUBDISTTRA + HUBDISTTR2 , data = pricemap, bw = lcrm3.bw,             
                      kernel = "bisquare", adaptive = TRUE, lambda.adjust = TRUE,cn.thresh = 30)
summary(lcrm3$SDF$Local_CN)

print(lcrm3)


test.CN <- function(model, data) {
    lcrmx.bw <- bw.gwr.lcr(model, data = data, kernel = "bisquare", adaptive = TRUE)
    print(model)
    print(lcrmx.bw)
    lcrmx <- gwr.lcr(model, data = data, bw = lcrmx.bw, kernel = "bisquare", adaptive = TRUE)
    print(summary(lcrmx$SDF$Local_CN))
    lcrmx$SDF$Local_CN}

data <- pricemap
model <- as.formula(CENA ~ AREA_SQ_M + HUBDISTBUS+ HUBDISTMET + HUBDISTTRA + HUBDISTTR2)
AllD <-test.CN(model, data)


model_no_train<- as.formula(CENA ~ AREA_SQ_M + HUBDISTBUS+ HUBDISTMET + HUBDISTTR2)
AllD_no_train <-test.CN(model_no_train, data)

model_no_tram<- as.formula(CENA ~ AREA_SQ_M + HUBDISTBUS+ HUBDISTMET + HUBDISTTRA)
AllD_no_tram <-test.CN(model_no_tram, data)


model_no_bus<- as.formula(CENA ~ AREA_SQ_M +  HUBDISTMET + HUBDISTTRA + HUBDISTTR2)
AllD_no_bus <-test.CN(model_no_bus, data)

model_no_metro<- as.formula(CENA ~ AREA_SQ_M + HUBDISTBUS+ HUBDISTTRA + HUBDISTTR2)
AllD_no_metro <-test.CN(model_no_metro, data)
summary(AllD_no_metro)

summary(AllD)
summary(AllD_no_train)
summary(AllD_no_tram)
summary(AllD_no_bus)
summary(AllD_no_metro)

#The model without train has the lowest condition number
lcrm_no_train.bw <- bw.gwr.lcr(CENA ~ AREA_SQ_M +  HUBDISTMET + HUBDISTBUS +
                         HUBDISTTR2 , data = pricemap, kernel = "bisquare", adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30)
lcrm_no_train.bw

lcrm_no_train <- gwr.lcr(CENA ~ AREA_SQ_M + HUBDISTMET + HUBDISTBUS + HUBDISTTR2 , data = pricemap, bw = lcrm_no_train.bw,             
                 kernel = "bisquare", adaptive = TRUE, lambda.adjust = TRUE,cn.thresh = 30)


print(lcrm_no_train)

lcrm_no_train_bus.bw <- bw.gwr.lcr(CENA ~ AREA_SQ_M +  HUBDISTMET + 
                                     HUBDISTTR2 , data = pricemap, kernel = "bisquare", adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30)

lcrm_no_train_bus <- gwr.lcr(CENA ~ AREA_SQ_M + HUBDISTMET + HUBDISTTR2 , data = pricemap, bw = lcrm_no_train_bus.bw,             
                         kernel = "bisquare", adaptive = TRUE, lambda.adjust = TRUE,cn.thresh = 30)
print(lcrm_no_train_bus)
#basic gwr
basic_no_train_bus.bw<- bw.gwr(CENA ~ AREA_SQ_M +  HUBDISTMET + 
                                  HUBDISTTR2 , data = pricemap, kernel = "bisquare", adaptive = TRUE)
basic_no_train_bus<-gwr.basic(CENA ~ AREA_SQ_M +  HUBDISTMET +  
                              HUBDISTTR2, data = pricemap,
                            bw = basic_no_train_bus.bw,kernel = "bisquare", adaptive = TRUE, F123.test = FALSE)
print(basic_no_train_bus)

#robust gwr
robust_no_train_bus.bw<- bw.gwr(CENA ~ AREA_SQ_M +  HUBDISTMET + 
                                 HUBDISTTR2 , data = pricemap, kernel = "bisquare", adaptive = TRUE)
robust_no_train_bus<-gwr.robust(CENA ~ AREA_SQ_M +  HUBDISTMET +  
                                HUBDISTTR2, data = pricemap,
                              bw = basic_no_train_bus.bw,kernel = "bisquare", adaptive = TRUE, F123.test = FALSE)

print(robust_no_train_bus)
#not good

#log transformation
log_basic_no_train_bus.bw<- bw.gwr(log(CENA) ~ AREA_SQ_M +  HUBDISTMET + 
                                 HUBDISTTR2 , data = pricemap, kernel = "bisquare", adaptive = TRUE)

log_basic_no_train_bus<-gwr.basic(log(CENA) ~ AREA_SQ_M +  HUBDISTMET +  
                                HUBDISTTR2, data = pricemap,
                              bw = log_basic_no_train_bus.bw,kernel = "bisquare", adaptive = TRUE, F123.test = FALSE)

print(log_basic_no_train_bus)
# not good
(summary(pricemap[,1:8]))


#map R2

results <-as.data.frame(basic_no_train_bus$SDF)
names(results)

summary(results$Local_R2)
gwr.map <- cbind(data, as.matrix(results))
library("tmap")
qtm(gwr.map, fill = "Local_R2")


mapRes <- tm_shape(gwr.map) + tm_fill("residual", n = 12, style = "jenks", palette = "Reds",
                                     title = "Residuals") +
  tm_layout(title= "Residuals sp. distr.", legend.text.size = 0.6, legend.title.size = 0.6, legend.position = c("RIGHT", "TOP"), frame =FALSE)

mapR2 <- tm_shape(gwr.map) + tm_fill("Local_R2", n = 12, style = "jenks", palette = "Reds",
                                    title = "Local R-square") +
  tm_layout(title= "Local R2 sp. distr.", legend.text.size = 0.6, legend.title.size = 0.6, legend.position = c("RIGHT", "TOP"), frame =FALSE)

#mapping coefficients

map1 <- tm_shape(gwr.map) + tm_fill("AREA_SQ_M", n = 12, style = "kmeans", palette = "Reds",
                                    title = "Area coefficient") +
  tm_layout(title= "Area coeff. sp. distr.", legend.text.size = 0.6, legend.title.size = 0.6, legend.position = c("RIGHT", "TOP"), frame =FALSE)

map2 <- tm_shape(gwr.map) + tm_fill("HUBDISTMET", n = 12, style = "kmeans", palette = "Reds",
                                    title = "Metro coefficient") +
  tm_layout(title= "Metro stop distance coeff. sp. distr.", legend.text.size = 0.6, legend.title.size = 0.6, legend.position = c("RIGHT", "TOP"), frame =FALSE)

map3 <- tm_shape(gwr.map) + tm_fill("HUBDISTTR2",  style = "fixed", breaks = c(-75, -50, -25, 0,25,50),palette = "Reds",
                                    title = "Tram coefficient") +
  tm_layout(title= "Tram stop distance coeff. sp. distr.", legend.text.size = 0.6, legend.title.size = 0.6, legend.position = c("RIGHT", "TOP"), frame =FALSE)


mypalette.10 <- brewer.pal(10, "Spectral")
spplot(basic_no_train_bus$SDF, "HUBDISTMET", key.space = "right",
          
            col.regions = mypalette.10, at = c(-45,-25,-10,-2,-1,0,1,2,10,30),
         
            main = "GWR coefficient estimates for Metro distance",
          
            sp.layout = map.layout)

summary(results$HUBDISTMET)
summary(results$HUBDISTTR2)


