clc
clear all
close all

syms s
A = [0 1;-6 -5]
B = [0;1]
C = [1 0]
D = [0]
sys=ss(A,B,C,D)
x_o = [0;1];

[num denum] = ss2tf(A,B,C,D) 
TF = C*inv(s.*eye(2)-A)*B + D


