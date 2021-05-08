c[n] = [1 .2 .4]

s[n] = 2*randi(2)-3
d[n] = sigma * randn()

s_hat[n] = ((s[n] * c[n]) + d[n]) * h[n]

Find R_sr[n], R_rr[n]

R_sr[n] = c[n] * R_ss[n]

r[n] = ((s[n] * c[n]) + d[n])

R_sr[n] = E[s(n)r(n-n0)] = E[s(n)(s[n-n0] * c[n-n0] + d[n-n0])]
= E[s(n)s(n-n0) * c[n-n0]] + E[s(n)d(n-n0)]
= E[S(n)s(n-n0)] * c[n]
= R_ss[n] * c[n]