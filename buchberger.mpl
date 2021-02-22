Buchberger:= proc(polys, variables)
    local i, j, f_i, f_j, g_i, g_j, G_, G_prime;
    G:={op(polys)}:
    while true do
        G_prime:=[op(G)];
        for i from 1 to nops(G_prime) do
            f_i := G_prime[i]:
            g_i := Groebner[LeadingMonomial](f_i, variables): 
            for j from i+1 to nops(G_prime) do
                f_j := G_prime[j]:
                g_j := Groebner[LeadingMonomial](f_j, variables):
                spoly_ := (lcm(g_i, g_j)/g_i)*f_i - (lcm(g_i, g_j)/g_j)*f_j:
                r := Groebner[Reduce](spoly_, G_prime, variables):
                if r<>0 then
                    G := G union {r}:
                end if:
            end do:
        end do:
        if G_prime=[op(G)] then
            break
        end if:
    end do:
    return G_prime:
end proc:


f:=[x^2+y^2-1, x*y-z, z^2-x^2-1, y^2+z^2-2];
variables := plex(z, y, x):

gb := Buchberger(f, variables):

printf("%a\n", gb):
# printf("%a\n", solve(gb, [z, y, x])):
f:=[x^2-x, x*y, y^2-y];
variables := plex(x, y):

gb := Buchberger(f, variables):

printf("%a\n", gb):
# gb := Groebner[Basis](f, variables):
# printf("%a\n", gb):
# printf("%a\n", solve(gb, [z, y, x])):