function u_snew = compute_u_snew(u_curr, u_prev, gamma)

u_snew = (gamma + 1)*u_curr - gamma*u_prev;

end