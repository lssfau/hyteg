# E := (g*(Rp^(-k)))/(nu*((k+1)*(k+2) - l*(l+1))*((k+3)*(k+4) - l*(l+1)));
Pp := (Ap*(r^l) + Bp*(r^(-l-1)) + Cp*(r^(l+2)) + Dp*(r^(-l+1)))*Y;
Pm := (Am*(r^l) + Bm*(r^(-l-1)) + Cm*(r^(l+2)) + Dm*(r^(-l+1)))*Y;

dPpdr := diff(Pp, r);
dPp2dr2 := diff(dPpdr, r);
dPp3dr3 := diff(dPp2dr2, r);

dPmdr := diff(Pm, r);
dPm2dr2 := diff(dPmdr, r);
dPm3dr3 := diff(dPm2dr2, r);

# # Freeslip
# eq||1 := subs(r = Rp, P)=0;
# eq||2 := subs(r = Rm, P)=0;
# eq||3 := subs(r = Rp, dP2dr2)=0;
# eq||4 := subs(r = Rm, dP2dr2)=0;

# # Noslip
# eq||1 := subs(r = R_plus, P)=0;
# eq||2 := subs(r = R_minus, P)=0;
# eq||3 := subs(r = R_plus, dPdr)=0;
# eq||4 := subs(r = R_minus, dPdr)=0;

# Both
eq||1 := subs(r = Rp, Pp)=0;
eq||2 := subs(r = Rm, Pm)=0;
eq||3 := subs(r = Rp, dPpdr)=0;
# eq||3 := subs(r = Rp, dPp2dr2)=0;
eq||4 := subs(r = Rm, dPm2dr2)=0;

# Interface conditions
eq||5 := subs(r = Rd, Pp - Pm)=0;
eq||6 := subs(r = Rd, dPpdr - dPmdr)=0;
eq||7 := subs(r = Rd, dPp2dr2 - dPm2dr2)=0;
eq||8 := subs(r = Rd, dPp3dr3 - dPm3dr3 - (g * Y/(nu * Rd)))=0;

ABCD := solve({eq||1, eq||2, eq||3, eq||4, eq||5, eq||6, eq||7, eq||8}, [Ap, Bp, Cp, Dp, Am, Bm, Cm, Dm]);
