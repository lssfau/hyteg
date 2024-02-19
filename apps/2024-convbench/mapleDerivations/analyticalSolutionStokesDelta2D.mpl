
psip := (Ap*(r^n) + Bp*(r^(-n)) + Cp*(r^(n+2)) + Dp*(r^(-n+2)))*sin(n*phi);
psim := (Am*(r^n) + Bm*(r^(-n)) + Cm*(r^(n+2)) + Dm*(r^(-n+2)))*sin(n*phi);

pp := ((-4*nu*Cp*(n+1))*r^n + (-4*nu*Dp*(n-1))*(r^(-n)) + F*(r^(k + 1)))*cos(n*phi);
pm := ((-4*nu*Cm*(n+1))*r^n + (-4*nu*Dm*(n-1))*(r^(-n)) + F*(r^(k + 1)))*cos(n*phi);

taup := nu*(diff(diff(psip, r), r) - ((1/r)*diff(psip, r)) - ((1/(r^2))*diff(diff(psip, phi), phi)));
taum := nu*(diff(diff(psim, r), r) - ((1/r)*diff(psim, r)) - ((1/(r^2))*diff(diff(psim, phi), phi)));

eq||1 := subs(r=Rp, diff(psip, phi))=0;
eq||2 := subs(r=Rm, diff(psim, phi))=0;

# eq||3 := subs(r=Rp, diff(psip, r))=0;
# eq||4 := subs(r=Rm, diff(psim, r))=0;

eq||3 := subs(r=Rp, taup)=0;
eq||4 := subs(r=Rm, taum)=0;

##### Interface conditions #####
dpsip := diff(psip, r);
dpsim := diff(psim, r);

d2psip := diff(dpsip, r);
d2psim := diff(dpsim, r);

d3psip := diff(d2psip, r);
d3psim := diff(d2psim, r);

eq||5 := subs(r=Rd, psip - psim)=0;
eq||6 := subs(r=Rd, dpsip - dpsim)=0;
eq||7 := subs(r=Rd, d2psip - d2psim)=0;
eq||8 := (subs(r=Rd, d3psip - d3psim) - (g*n*sin(n*phi)/(nu*Rd)))=0;
# eq||8 := subs(r=Rd, pp - pm)=0;

ABCD := solve({eq||1, eq||2, eq||3, eq||4, eq||5, eq||6, eq||7, eq||8}, [Ap, Bp, Cp, Dp, Am, Bm, Cm, Dm]);

# with(CodeGeneration);

# C(ABCD[1][1]);
# C(ABCD[1][2]);
# C(ABCD[1][3]);
# C(ABCD[1][4]);
