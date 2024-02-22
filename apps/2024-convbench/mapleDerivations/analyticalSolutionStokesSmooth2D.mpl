
psi := (A*(r^n) + B*(r^(-n)) + C*(r^(n+2)) + D*(r^(-n+2)) + E*(r^(k+3)))*sin(n*phi);
tau := nu*(diff(diff(psi, r), r) - ((1/r)*diff(psi, r)) - ((1/(r^2))*diff(diff(psi, phi), phi)));

eq||1 := subs(r=Rp, diff(psi, phi))=0;
eq||2 := subs(r=Rm, diff(psi, phi))=0;

eq||3 := subs(r=Rp, diff(psi, r))=0;
# eq||4 := subs(r=Rm, diff(psi, r))=0;

# eq||3 := subs(r=Rp, tau)=0;
eq||4 := subs(r=Rm, tau)=0;

ABCD := solve({eq||1, eq||2, eq||3, eq||4}, [A, B, C, D]);

with(CodeGeneration);

# C(ABCD[1][1]);
# C(ABCD[1][2]);
# C(ABCD[1][3]);
# C(ABCD[1][4]);
