### replicate Herwartz,Wang 2024:630, Ch.4 ###
\donttest{
# select minimal or full example #
is_min = TRUE
n.boot = ifelse(is_min, 5, 1000)
idx_i  = if(is_min){ 1 }else{ 1:14 }

# load and prepare data #
data("EURO")
names_i = names(EURO[idx_i+1])  # country names (#1 is EA-wide aggregated data)
names_s = paste0("epsilon[ ", c(1:2, "m", "f"), " ]")  # shock names
idx_k   = 1:4   # endogenous variables in individual data matrices
idx_t   = 1:(nrow(EURO[[1]])-1)  # periods from 2001Q1 to 2019Q4 
trend2  = idx_t^2

# panel SVARX model, Ch.4.1 #
L.data = lapply(EURO[idx_i+1], FUN=function(x) x[idx_t, idx_k])
L.exog = lapply(EURO[idx_i+1], FUN=function(x) cbind(trend2, x[idx_t, 5:10]))
R.lags = c(1,2,1,2,2,2,2,2,1,2,2,2,2,1)[idx_i]; names(R.lags) = names_i
R.pvar = pvarx.VAR(L.data, lags=R.lags, type="both", D=L.exog)
R.pid  = pid.dc(R.pvar, combine="pool")
print(R.pid)  # suggests e3 and e4 to be MP and financial shocks, respectively.
colnames(R.pid$B) = names_s  # accordant labeling

# EA-wide SVARX model, Ch.4.2 #
R.data = EURO[[1]][idx_t, idx_k]
R.exog = cbind(trend2, EURO[[1]][idx_t, 5:6])
R.varx = pvarx.VAR(list(EA=R.data), lags=2, type="both", D=list(EA=R.exog))
R.id   = pid.dc(R.varx, combine="indiv")$L.varx$EA
colnames(R.id$B) = names_s  # labeling

# comparison of IRF without confidence bands, Ch.4.3.1 #
data("EU_w")  # GDP weights with the same ordering names_i as L.varx in R.pid
R.norm = function(B) B / matrix(diag(B), nrow(B), ncol(B), byrow=TRUE) * 25
R.irf  = as.pplot(
  EA=plot(irf(R.id,  normf=R.norm), selection=list(idx_k, 3:4)),
  MG=plot(irf(R.pid, normf=R.norm, w=EU_w), selection=list(idx_k, 3:4)),
  color_g=c("#3B4992FF", "#008B45FF"), shape_g=16:17, n.rows=length(idx_k))
plot(R.irf)

# comparison of IRF with confidence bands, Ch.4.3.1 #
R.boot_EA = sboot.mb(R.id, b.length=8, n.boot=n.boot, n.cores=2, normf=R.norm)
R.boot_MG = sboot.pmb(R.pid, b.dim=c(8, FALSE), n.boot=n.boot, n.cores=2, 
                      normf=R.norm, w=EU_w)
R.irf = as.pplot(
  EA=plot(R.boot_EA, lowerq=0.16, upperq=0.84, selection=list(idx_k, 3:4)),
  MG=plot(R.boot_MG, lowerq=0.16, upperq=0.84, selection=list(idx_k, 3:4)),
  color_g=c("#3B4992FF", "#008B45FF"), shape_g=16:17, n.rows=length(idx_k))
plot(R.irf)
}

