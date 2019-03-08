# Algorithm-PSIA

This is for the paper: Sun, Yahui, et al. "A physarum-inspired prize-collecting steiner tree approach to identify subnetworks for drug repositioning." BMC systems biology 10.5 (2016): 128.

This repository contains the MATLAB coding of the PSIA algorithm, which is used to identify subnetworks in drug similarity network.

There are five MATLAB files in total: PSIA.m is the main coding for the PSIA algorithm; D_01_b.mat is the data of the input node-weighted graph; Function_OutputDegree.m is the function used to calculate the degrees of vertices; Function_prim.m is the function used to delete edges which are not in the same connected component with the terminals; Function_SMTpost.m is the function used to find the Minimum Spanning Tree of the input graph.

To run PSIA, put all the five MATLAB files in the same document, and run PSIA.m.

For any issue, please feel free to email me, syhhit@gmail.com   -Yahui

