# E := (g*(Rp^(-k)))/(nu*((k+1)*(k+2) - l*(l+1))*((k+3)*(k+4) - l*(l+1)));
P := (A*(r^l) + B*(r^(-l-1)) + C*(r^(l+2)) + D*(r^(-l+1)) + E*(r^(k+3)))*Y;
dPdr := diff(P, r);
dP2dr2 := diff(dPdr, r);

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
eq||1 := subs(r = Rp, P)=0;
eq||2 := subs(r = Rm, P)=0;
eq||3 := subs(r = Rp, dPdr)=0;
eq||4 := subs(r = Rm, dP2dr2)=0;

ABCD := solve({eq||1, eq||2, eq||3, eq||4}, [A, B, C, D]);
