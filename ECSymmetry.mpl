with(LinearAlgebra):
with(combinat):
with(Iterator):
with(ListTools):

esvp := proc(f,v)   
# Extended Sign Variation of Polynomial
# In : f, a polynomial in v with real coefficients
#      v, a variable
# Out: the sign variataion count on the coefficients of f
  local d,L,Lb, Cb, i,c;
  d  := degree(f,v);
  L  := [coeff(f,v,d-i)$i=0..d];
  c  := esvl(L);
  return c;
end:

esvl := proc(L)  
# Extended Sign Variation of List
# In : L, a list of real numbers
# Out: extended sign variation count on L
  local x,cls;
  cls := invclosure(L);
  return add(svl(x), x in cls) / nops(cls);
end:

invclosure := proc(L)
# Sign sequence inverse closure
# In : L, a list of signs
# Out: L', list of elements obtained from L where each 0 is replaced with + and -
  local inds,es,svcs,all_seqs,e,Le,i;
  inds := [SearchAll(0, L)];
  es := cpl([[-1,1]$nops(inds)]);
  svcs := 0;
  all_seqs := [];
  for e in es do
    Le := L;
    for i from 1 to nops(inds) do
      Le[inds[i]] := e[i];
    od:
    all_seqs := [op(all_seqs), Le];
  od:
  return all_seqs;
end:
 
svl := proc(L)    
# Sign Variation of List
# In : L, a list of real numbers
# Out: the sign variation count on L
  local i,c,v,C;
  C := [seq(L[i]*L[i+1],i=1..nops(L)-1)];
  C := remove(v->evalf(v>0), C);
  c  := nops(C);
  return c;   
end:

cpl := proc(Ls)    
# Cartesian Product of Lists
# In : Ls, a list of lists
# Out: The cartesian product of Ls
  local P,p,C;
  P := cartprod(Ls);
  C := [];
  while not P[finished] do
    C := [op(C), P[nextvalue]()];
  end do;
  
  return C;
end:

ESP := proc(k, X)
# Elementary Symmetric Polynomial
# In  : k, integer between 1 and n
#       X, list of variables
# Out : the kth symmetric polynomial in X
  local s,x;
  
  return add(mul(x, x in s), s in choose(X, k)):
end:

FTSP := proc(h, X, ov)
# Fundamental Theorem of Symmetric Polynomials
# In  : h, symmetric polynomial in X
#     : X, list of variables
#     : ov, output variable
# Out : p, polynomial in ov[1],...,ov[n] such that h = p(s(1,n), ..., s(n, n)), 
#          where s(k, n) are the elementary symmetric polynomials on X
  local n, hs, ps, g, i, d, w, p;

  n := nops(X);
  
  # base case: n = 1 or totaldegree(h) = 0
  if n = 1 or degree(h, {op(X)}) <= 0 then
    return eval(h, [seq(X[i] = ov[i], i=1..n)]);
  fi;
  
  hs := eval(h, X[n] = 0);
  ps := FTSP(hs, X[1..n-1], ov);
  g := eval(ps, [seq(ov[i] = ESP(i, X), i=1..n-1)]);
  d := expand((h - g)/ESP(n, X));
  w := FTSP(d, X, ov);
  p := w * ov[n] + ps;
  return expand(p);
end:

TM := proc(m)
# T Matrix
# In : m, positive integer
# Out: T, T matrix
  local Tm,i,j,t,k,r,iss,is;
  Tm := Matrix(m,m);
  for k from 1 to m do
    iss := choose(m,k);
    for r from 1 to m do
      t := 0;
      for is in iss do
        if irem(nops(select(e->e<=r,is)),2)=1 then
          t := t + 1; 
        fi;
      od;
      Tm[k,r] := t;
    od;
  od;
  return Tm;
end:

Csym := proc(m)
# Combinatorial part, symmetry method
# In : m, positive integer
# Out: Csym
  return TM(m)^(-1);
end:

hv := proc(m,n,r)
# h Polynomial
# In : m, n, r, positive integers
# Out: hr, the r-th h polynomial
  local hr,iss,is,i,j,w;
  hr := 1;
  iss := choose(m,r);
  for j from 1 to n do
    for is in iss do
      w := mul(alpha[i]-beta[j],i=is);
      hr := hr*(x+w);
    od;
  od;
  return hr;
end:

Asym := proc(F,G)
# Algebraic part, symmetric polynomial method
# In : F, G, symbolic real symmetric matrices
# Out: Asym, a vector of polynomials
  local m,n,f,g,As,r,p,h,i;
  m := RowDimension(F);
  n := RowDimension(G);
  f := CharacteristicPolynomial(F, x);
  g := CharacteristicPolynomial(G, x);

  As := [];
  for r from 1 to m do
    h := hv(m, n, r);
    p := FTSP(h, [seq(beta[i],   i=1..n)], y);
    p := FTSP(p, [seq(alpha[i],  i=1..m)], z);
    p := eval(p, [seq(y[i]=coeff(g,x,n-i) * (-1)^i, i=1..n)]);
    p := eval(p, [seq(z[i]=coeff(f,x,m-i) * (-1)^i, i=1..m)]);
    As := [op(As), collect(p,x)];
  od:

  return Vector(As);
end:

EC := proc(F,G)    
# Eigenvalue Configuration, symmetric polynomial method
# In : F, G, numeric real symmetric matrices
# Out: Eigenvalue configuration of F and G   
  local m,Tm,Dsc,vd,d,c;
  m   := RowDimension(F);
  Tm  := TM(m);
  Dsc := Dsym(F,G);
  vd  := Vector([seq(esvp(d,x), d in Dsc)]);
  c   := Tm^(-1) . vd;
  return c;
end:

ECSym := proc(F,G)
# Conditions for EC, symmetric polynomial method
# In : F, G, symbolic real symmetric matrices
# Out: Csym, Asym such that c = EC(F,G) iff c = Csym . v(Asym)
  return Csym(RowDimension(F)), Asym(F,G);
end:
