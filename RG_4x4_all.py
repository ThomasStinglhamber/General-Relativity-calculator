#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 14:56:53 2020

@author: Thomas Stinglhamber (De La Physique)
"""


import sympy as sp
import numpy as np
from pprint import pprint

""" definition constantes """

t = sp.Symbol('t')
r = sp.Symbol('r')
o= sp.Symbol('o') # theta
p= sp.Symbol('p') # phi
G= sp.Symbol('G')
M= sp.Symbol('M')
u= sp.Symbol('u')
v= sp.Symbol('v')
K= sp.Symbol('K')
p_1 = sp.Symbol('p_1')
p_2 = sp.Symbol('p_2')
p_3 = sp.Symbol('p_3')
p_4 = sp.Symbol('p_4')
C = sp.Symbol('C')

""" Definition metrique """ 

g = [[((1-2*G*M)/r) ,0,0,0],
     [0,(1/(1-2*G*M/r)),0,0],
     [0,0,r**2,0],
     [0,0,0,r**2 * sp.sin(o)**2]]


""" variable de la métrique """


k=[t,r,o,p]

""" Chemin pour créer le fichier .txt avec le code LaTeX """

Path ="/Users/thomasstinglhamber/Desktop/RG_ALL/LaTeX_RG_4x4.txt"# changer pour mettre votre chemin


ToLaTeX = True # False pour pas générer le fichier.txt


#--------------------- Plus rien a modifier à partir d'ici -----------------------
f = open(Path, "w")

g2=sp.Matrix(([g[0][0],g[1][0],g[2][0],g[3][0]],
                      [g[0][1],g[1][1],g[2][1],g[3][1]],
                      [g[0][2],g[1][2],g[2][2],g[3][2]],
                      [g[3][0],g[3][1],g[3][2],g[3][3]]))
g_inv= g2.inv()
        

print("\u0332".join("Metrique : "))
pprint(g)
print()

    
def Christo4x4():
    Christo = [[0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],]
    Christo_t = [[0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],]
    Christo_r = [[0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],]
    Christo_o = [[0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],]
    Christo_p = [[0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],]
    Christofell=[]
        
    """ definition des symboles de Christofell """
    
    
    print("\u0332".join("Symbole de Christofell: "))
    for m in range(4):
        for l in range(4):
            for i in range(4):
                for j in range(4):
                    Christo[i][j] = 1/2 * g_inv[l,m]*(sp.diff(g[l][i],k[j])+sp.diff(g[j][l],k[i])-sp.diff(g[i][j],k[l]))
                    if Christo[i][j] !=0:
                        print("C|",k[m],k[i],k[j],":",sp.simplify(Christo[i][j]) )   
                        if m ==0:
                            Christo_t[i][j] =Christo[i][j]
                        if m ==1 :
                            Christo_r[i][j] =Christo[i][j]
                        if m ==2 :
                            Christo_o[i][j] =Christo[i][j]
                        if m ==3 :
                            Christo_p[i][j] =Christo[i][j]
                        
    Christofell =[Christo_t,Christo_r,Christo_o,Christo_p]
    print()
    for i in range(4):
        print("Christofell selon",k[i],":")
        pprint(Christofell[i])
        print()
        
    return Christofell


def Riemann(Christofell):
    
    """ defintion de matrice utile apres pour indexation"""
    
    Ri = [[0,0,0,0],
          [0,0,0,0],
          [0,0,0,0],
          [0,0,0,0],]
    Rl0=[[0,0,0,0],
          [0,0,0,0],
          [0,0,0,0],
          [0,0,0,0],]
    R_l0n0=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l0n1=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l0n2=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l0n3=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    Rl1=[[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],]
    R_l1n0=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l1n1=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l1n2=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l1n3=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l2n0=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l2n1=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l2n2=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l2n3=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    Rl2=[[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],]
    R_l3n0=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l3n1=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l3n2=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    R_l3n3=[[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],]
    Rl3=[[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],]
    R =[]
    
    """ définition du tenseur de Riemann """
    
    print("\u0332".join("Tenseur de Riemann :"))
    for l in range(4):
        for i in range(4):
            for j in range(4):
                for n in range(4):
                    for m in range(4):
                        Ri[j][n] = ((sp.diff(Christofell[l][i][j],k[n]))-(sp.diff(Christofell[l][i][n],k[j])) +((Christofell[m][i][j])*(Christofell[l][m][n]))-((Christofell[m][i][n])*(Christofell[l][m][j])))
                    if m==3:
                        if Ri[j][n] !=0:
                            print("R|",k[l],k[i],k[j],k[n],":",sp.simplify(Ri[j][n]))
                            if l==0:
                                if i ==0:
                                    R_l0n0[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l0n1[j][n]= Ri[j][n]
                                if i ==2:
                                    R_l0n2[j][n]= Ri[j][n]
                                if i ==2:
                                    R_l0n3[j][n]= Ri[j][n]
                            if l==1:
                                if i ==0:
                                    R_l1n0[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l1n1[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l1n2[j][n]= Ri[j][n]
                                if i ==2:
                                    R_l1n3[j][n]= Ri[j][n]
                            if l==2:
                                if i ==0:
                                    R_l2n0[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l2n1[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l2n2[j][n]= Ri[j][n]
                                if i ==2:
                                    R_l2n3[j][n]= Ri[j][n]
                            if l==3:
                                if i ==0:
                                    R_l3n0[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l3n1[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l3n2[j][n]= Ri[j][n]
                                if i ==2:
                                    R_l3n3[j][n]= Ri[j][n]
    
    Rl0 =[R_l0n0,R_l0n1,R_l0n2,R_l0n3]
    Rl1 =[R_l1n0,R_l1n1,R_l1n2,R_l1n3]
    Rl2 =[R_l2n0,R_l2n1,R_l2n2,R_l2n3]
    Rl3 =[R_l3n0,R_l3n1,R_l3n2,R_l3n3]
    R =[Rl0,Rl1,Rl2,Rl3]
    print()
    for i in range(4):  
        print("Riemann selon",k[i],":")
        pprint(R[i])
        print()                 
    return R 

def Ricci(R):

    ricci1=[[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],]
    ricci=[[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],]

    """ Ricci """
    print("\u0332".join("Tenseur de Ricci :"))
    for l in range(4):
        for i in range(4):
            for j in range(4):
                ricci1[i][j] = R[l][i][j][l]
                if ricci1[i][j]!=0:
                    ricci[i][j]=ricci1[i][j]
                        
    print(ricci)
    print()
    return ricci

def Courbure_scalaire(ricci):
    scalaire=[[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],]
    
    """ courbure  scalaire """
    
    print("\u0332".join("Courbure scalaire :"))
    for i in range(4):
        for j in range(4):
            scalaire[i][j]= g_inv[i,j]* ricci[i][j]
            
    somme = scalaire[0][0]+scalaire[1][0]+scalaire[2][0]+scalaire[3][0]+scalaire[0][1]+scalaire[1][1]+scalaire[2][1]+scalaire[3][1]+scalaire[0][2]+scalaire[1][2]+scalaire[2][2]+scalaire[3][2]+scalaire[0][3]+scalaire[1][3]+scalaire[2][3]+scalaire[3][3]
    print(somme)
    return somme


"""--------------- Partie pour générer le code LaTeX ----------------- """


def PythonToLaTeX(Christofell,R,ricci,somme):
        
    def pmatrix(a):
        
        lines = str(a).replace('[', '').replace(']', '').replace(' - ','-').replace(' + ','+').splitlines()
        rv = [r' $$\begin{pmatrix}']
        rv += ['  ' + ' ,& '.join(l.split()) + r'\\' for l in lines]
        rv +=  [r'\end{pmatrix}$$']   
        
        return '\n'.join(rv)
    
    def Tensor_matrix(a):
        rv = [r'$$ \begin{pmatrix}']
        for i in range(len(a)):
            rv += [r' \begin{pmatrix}']
            for j in range(len(a[i])):
                lines = str(a[i][j]).replace('[', '').replace(']', '').replace(' - ','-').replace(' + ','+').splitlines()
                rv += ['  ' + ', & '.join(l.split()) + r'\\' for l in lines]
            rv +=  [r'\end{pmatrix}\\'] 
        rv +=  [r'\end{pmatrix}$$'] 
        return '\n'.join(rv)
    
    print('\setlength\parindent{0pt}',file=f)
    print('\section{RG4x4}',file=f)
    print('\subsection{Métrique}',file=f)
    print(pmatrix(np.array(g)),'\n',file=f)
    print('\subsection{Christofell}',file=f)
    for m in range(4):
        for i in range(4):
            for j in range(4):
                if Christofell[m][i][j] != 0:
                    print('$\Gamma^{',k[m],'}_{',k[i],k[j],'}',"=",Christofell[m][i][j],r'$\\', file=f)
                else :pass
        print('\n',file=f)
            
        
# =============================================================================
#     for i in range(4):
#         print('Symbole de Christofell selon $',k[i],'$ :\n',pmatrix(np.array(Christofell[i])) + '\n', file=f)
#     
# =============================================================================
    print('\subsection{Riemann}',file=f)
    for l in range(4):
        for i in range(4):
            for j in range(4):
                for n in range(4):
                    if R[l][i][j][n] !=0:
                        print('$R^{',k[l],'}_{',k[i],k[j],k[n],'}',"=",R[l][i][j][n],r'$\\', file=f)
        print('\n',file=f)
        
        
# =============================================================================
#     for i in range(4):
#         print('Tenseur de Riemann selon $',k[i],'$ :\n',Tensor_matrix(np.array(R[i])) + '\n', file=f)
#         
# =============================================================================
    print('\subsection{Ricci}',file=f)
    for l in range(4):
        for i in range(4):
            for j in range(4):
                if R[l][i][j][l] !=0:
                        print('$R\mathrm{icci}_{',k[i],k[j],'}',"=",R[l][i][j][l],r'$\\', file=f)
        print('\n',file=f)
    
# =============================================================================
#     print('Tenseur de Ricci ','\n',pmatrix(np.array(ricci)) + '\n', file=f)
#     
# =============================================================================
    print('\paragraph{Courbure scalaire} ','\n','$$',somme , '$$',file=f)
    
    f.close()


Christofell = Christo4x4()
R =Riemann(Christofell)
ricci =Ricci(R)
somme =Courbure_scalaire(ricci) 

if ToLaTeX == True:
    PythonToLaTeX(Christofell,R,ricci,somme) 
    f1 = open(Path, 'r')
    input_data = f1.read()
    f1.close()
    
    theta_1 =' \sin(\\theta)'
    theta_4 =' \cos(\\theta)'
    theta_2 =' \\theta '
    phi_2 =' \sin(\phi)'
    phi_3 =' \cos(\phi)'
    input_data = input_data.replace('o ', theta_2.replace('\t', '')).replace('p ' , ' \phi ').replace('**', '^').replace('cos(o)', theta_4.replace('\t', '')).replace('sin(o)', theta_1.replace('\t', '')).replace('sin(p)', phi_2).replace('cos(p)', phi_3).replace('*' , '').replace('1.00000000000000' , '1').replace('1.0/' , '1/').replace('1.0+' , '1+').replace('1.0-' , '1-').replace('1.0' , '')                                          
    
    f2 = open(Path, 'w')
    f2.write(input_data)
    f2.close()
    print('\n','Fichier .txt généré')
else :pass

