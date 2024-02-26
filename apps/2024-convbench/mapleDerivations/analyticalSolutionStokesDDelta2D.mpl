
psip := (Ap*(r^n) + Bp*(r^(-n)) + Cp*(r^(n+2)) + Dp*(r^(-n+2)))*sin(n*phi);
psipm := (Apm*(r^n) + Bpm*(r^(-n)) + Cpm*(r^(n+2)) + Dpm*(r^(-n+2)))*sin(n*phi);
psim := (Am*(r^n) + Bm*(r^(-n)) + Cm*(r^(n+2)) + Dm*(r^(-n+2)))*sin(n*phi);

taup := nu*(diff(diff(psip, r), r) - ((1/r)*diff(psip, r)) - ((1/(r^2))*diff(diff(psip, phi), phi)));
taupm := nu*(diff(diff(psipm, r), r) - ((1/r)*diff(psipm, r)) - ((1/(r^2))*diff(diff(psipm, phi), phi)));
taum := nu*(diff(diff(psim, r), r) - ((1/r)*diff(psim, r)) - ((1/(r^2))*diff(diff(psim, phi), phi)));

eq||1 := subs(r=Rp, diff(psip, phi))=0;
eq||2 := subs(r=Rm, diff(psim, phi))=0;

# eq||3 := subs(r=Rp, diff(psip, r))=0;
# eq||4 := subs(r=Rm, diff(psim, r))=0;

eq||3 := subs(r=Rp, taup)=0;
eq||4 := subs(r=Rm, taum)=0;

##### Interface conditions #####
dpsip := diff(psip, r);
dpsipm := diff(psipm, r);
dpsim := diff(psim, r);

d2psip := diff(dpsip, r);
d2psipm := diff(dpsipm, r);
d2psim := diff(dpsim, r);

d3psip := diff(d2psip, r);
d3psipm := diff(d2psipm, r);
d3psim := diff(d2psim, r);

eq||5 := subs(r=Rdp, psip - psipm)=0;
eq||6 := subs(r=Rdp, dpsip - dpsipm)=0;
eq||7 := subs(r=Rdp, d2psip - d2psipm)=0;
eq||8 := (subs(r=Rdp, d3psip - d3psipm) - (g*n*sin(n*phi)/(nu*Rdp)))=0;

eq||9 := subs(r=Rdm, psipm - psim)=0;
eq||10 := subs(r=Rdm, dpsipm - dpsim)=0;
eq||11 := subs(r=Rdm, d2psipm - d2psim)=0;
eq||12 := (subs(r=Rdm, d3psipm - d3psim) - (g*n*sin(n*phi)/(nu*Rdm)))=0;

ABCD := solve({eq||1, eq||2, eq||3, eq||4, eq||5, eq||6, eq||7, eq||8, eq||9, eq||10, eq||11, eq||12}, [Ap, Bp, Cp, Dp, Apm, Bpm, Cpm, Dpm, Am, Bm, Cm, Dm]);

# with(CodeGeneration);

# C(ABCD[1][1]);
# C(ABCD[1][2]);
# C(ABCD[1][3]);
# C(ABCD[1][4]);
