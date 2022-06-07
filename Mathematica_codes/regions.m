(* This code implements a simple Monte Carlo analysis of feasibility
  regions as described in section S4  of "Understanding the
  variability of pairwise coexistence within multispecies systems" by
  Deng, Taylor and Saavedra.  Code by WT, January 2022 *)

(* 
  To run the code, input into a mathematica shell using
  "<< regions.m"
  Defining a matrix mx,  "statistics[mx, n]" generates a list of
  frequencies for the different subset feasibilities for the matrix mx
  with n separate random samples on the unit sphere.  Subset number
  indexes into "subsets[d]", which gives all subsets in d dimensions

Sample run:

% compute
Mathematica 11.1.1 Kernel for Linux x86 (64-bit)
Copyright 1988-2017 Wolfram Research, Inc.

In[1]:= << regions.m

In[2]:= mx ={{1, 0.2, 0.1},{0.2, 1, 0.3},{0.1, 0.3, 1}}

Out[2]= {{1, 0.2, 0.1}, {0.2, 1, 0.3}, {0.1, 0.3, 1}}

In[3]:= statistics[mx, 10000]

Out[3]= {0.1496, 0.1763, 0.1552, 0.1074, 0.1403, 0.0878, 0.0542}

In[4]:= statistics[mx, 10000]

Out[4]= {0.1552, 0.1662, 0.1573, 0.117, 0.1369, 0.0889, 0.0543}

In[5]:= subsets[3]

Out[5]= {{1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}}

So for example the size of the feasibility region for all 3 species is
roughly 0.54, while the feasibility of only species 1 is around 0.15.
Running with larger n gives higher accuracy.  Note that these
feasibility regions are normalized relative to the complete sphere
including negative values of the parameters.

 *)

(* The main component of the code implements the algorithm that *)
(* determines boundaries between regions with different feasible sets *)

(* selects nonzero elements of a vector  associated with subset slots *)


extract[d_, v_, slots_] :=
  Table[If[MemberQ[slots, i], v[[i]], 0],{i, 1, d}]

(* unit vector in direction i *)

unit[d_, i_] :=Table[If[j == i, 1, 0],{j, 1, d}]

(* computes determinant of relevant matrix (details in paper) *)

(* d = dimension
   m = interaction matrix
   s =(smaller) subset of directions
   i = index added to s for orientation 
   v = vector being tested*)

(* this will be positive if v is on the same side as i unit vector *)

side[d_, m_, s_, i_, v_] :=
  Det[Table[If[MemberQ[s, j], extract[d, m[[j]], Join[s,{i}]],
            If[j == i, v, unit[d, j]]],{j, 1, d}]]

(* same but also if i is in s, gives proper sign for s subset *)

sides[d_, m_, s_, i_, v_]:=
  If[MemberQ[s, i], side[d, m, Complement[s,{i}], i, v],
        -side[d, m, s, i, v]]

(* streamline code so don't have to keep computing, by computing
  linear coefficients for conditions *)

coefficients[d_, m_, s_, i_] :=
 coefficients[d, m, s, i]=
  Block[{ss=  sides[d, m, s, i,Table[x[k],{k, 1, d}]]},
    Table[Coefficient[ss, x[j]],{j, 1, d}]]

(* determine region by simple linear relations *)

region[d_, m_, s_, v_] :=
  Apply[And, 
    Table[(v.coefficients[d, m, s, i]> 0),{i, 1, d}]]

(* now compute statistics *)

(* gives all non-empty subsets in dimension d *)

subsets[d_] := Drop[Subsets[Table[i,{i, 1, d}]], 1]

(* checks all conditions for a given vector *)

conditions[d_, m_, v_] :=
 Map[If[region[d, m, #, v], 1, 0]&, subsets[d]]

(* finds the subset associated with the feasibility domain for matrix m *)
(* choosing each x from a normal distribution gives x uniformly
  distributed on the unit sphere. *)

(* working in positive orthant *)
random[d_, m_] :=
  conditions[d, m, Map[Abs, RandomReal[ NormalDistribution[0, 1],{d}]]]

(* without taking absolute value *)
random[d_, m_] :=
  conditions[d, m, RandomReal[ NormalDistribution[0, 1],{d}]]



(* doing statistics on multiple random samples *)

multiple[d_, m_, n_] := 
  Apply[Plus,Table[random[d, m],n]]

statistics[m_, n_] :=
 With[{results = multiple[Length[m], m, n]},
   N[results/n]]
