function f = lower_confidence_bound(mu,stdv,kappa)
 f = mu - kappa*stdv;
end