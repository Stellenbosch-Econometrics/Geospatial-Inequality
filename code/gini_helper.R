

gini_wiki <- function(x) {1 + 2/(length(x)-1) * (sum(seq_along(x)*sort(x)) / sum(x) - length(x))} 

gini_noss <- function(x) 2/length(x) * sum(seq_along(x)*sort(x)) / sum(x) - (length(x)+1)/length(x) # same as in Galimberti et al. (2020)

# S <- function(p, q) {q/2 * (2*p + q - 1)}
# all.equal(sum(30:100), S(30, 70+1))
# Sn <- function(n) n*(n+1)/2
# all.equal(sum(1:100), Sn(100))
# Skn <- function(k, n) (n*(n+1) - (k-1)*k)/2 
# all.equal(Skn(30, 100), sum(30:100))
# Skpq <- function(k, q) q/2 * (2*k + q + 1) + k
# all.equal(Skpq(30, 70), sum(30:100))
Skp1qm1 <- function(k, q) (q-1)/2 * (2*(k+1) + q) + k + 1
all.equal(Skp1qm1(30-1, 70+1), sum(30:100))

w_gini <- function(x, w, sscor = FALSE) {
  si = .Internal(qsort(x, TRUE))
  w = w[si$ix]
  x = si$x
  sw = sum(w)
  csw = cumsum(w)
  sx = Skp1qm1(c(0, csw[-length(csw)]), w) # Skpq(c(0, csw[-length(csw)])+1, w-1)
  if(sscor) return(1 + 2/(sw-1)*(sum(sx*x) / sum(x*w) - sw)) # TODO: see paper again, correction should be different 
  2/sw * sum(sx*x) / sum(x*w) - (sw+1)/sw
}
