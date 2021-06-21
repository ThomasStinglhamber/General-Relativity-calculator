#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 14:03:09 2020

@author: Thomas Stinglhamber (De La Physique)
"""

import sympy as sp
import numpy as np
from pprint import pprint

""" definition constantes """
r = sp.Symbol('r')
o = sp.Symbol('o') # theta
p = sp.Symbol('p') # phi

s = sp.Symbol('s')
z = sp.Symbol('z')

""" Definition metrique """ 
g = [[1 ,0,0],
     [0,r**2 ,0],
     [0,0,r**2 * sp.sin(o)**2]]


""" variable de la métrique """
k=[r,o,p]                  # on derive la metrique par rapport a ces indices


""" Chemin pour créer le fichier .txt avec le code LaTeX """

Path ="/Users/thomasstinglhamber/Desktop/RG_ALL/LaTeX_RG_3x3.txt"# changer pour mettre votre chemin


ToLaTeX = True # False pour pas générer le fichier.txt


#--------------------- Plus rien a modifier à partir d'ici -----------------------
f = open(Path, "w")

g2=sp.Matrix(([g[0][0],g[1][0],g[2][0]],
              [g[0][1],g[1][1],g[2][1]],
              [g[0][2],g[1][2],g[2][2]]))

g_inv= g2.inv()



print("\u0332".join("Metrique : "))
pprint(g)
print()

def Christo3x3():
    """ defintion de matrice utile apres pour indexation"""
    Christo = [[0,0,0],
               [0,0,0],
               [0,0,0],]
    Christo_1 = [[0,0,0],
                 [0,0,0],
                 [0,0,0],]
    Christo_2 = [[0,0,0],
                 [0,0,0],
                 [0,0,0],]
    Christo_3 = [[0,0,0],
                 [0,0,0],
                 [0,0,0],]
    
    
    """ definition des symboles de Christofell """
    
    
    print("\u0332".join("Symbole de Christofell: "))
    for m in range(3):
        for l in range(3):
            for i in range(3):
                for j in range(3):
                    Christo[i][j] = 1/2 * g_inv[l,m]*(sp.diff(g[l][i],k[j])+sp.diff(g[j][l],k[i])-sp.diff(g[i][j],k[l]))
                    if Christo[i][j] !=0:
                        print("C|",k[m],k[i],k[j],":",Christo[i][j])    
                        if m ==0:
                            Christo_1[i][j] =Christo[i][j]
                        if m ==1 :
                            Christo_2[i][j] =Christo[i][j]
                        if m ==2 :
                            Christo_3[i][j] =Christo[i][j]
    Christofell=[]
    Christofell =[Christo_1,Christo_2,Christo_3]
    print()
    for i in range(3):
        print("Christofell selon",k[i],":")
        pprint(Christofell[i])
        print()
    
    return Christofell



def Riemann(Christofell):
        
    """ defintion de matrice utile apres pour indexation"""
    Ri = [[0,0,0],
          [0,0,0],
          [0,0,0],]
    Rl0=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R_l0n0=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R_l0n1=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R_l0n2=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    Rl1=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R_l1n0=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R_l1n1=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R_l1n2=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R_l2n0=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R_l2n1=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R_l2n2=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    Rl2=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    R =[]
    
    """ définition du tenseur de Riemann """
    
    print("\u0332".join("Tenseur de Riemann :"))
    for l in range(3):
        for i in range(3):
            for j in range(3):
                for n in range(3):
                    for m in range(3):
                        Ri[j][n] = ((sp.diff(Christofell[l][i][j],k[n]))-(sp.diff(Christofell[l][i][n],k[j])) +((Christofell[m][i][j])*(Christofell[l][m][n]))-((Christofell[m][i][n])*(Christofell[l][m][j])))
                    if m==2:
                        if Ri[j][n] !=0:
                            print("R|",k[l],k[i],k[j],k[n],":",(Ri[j][n]))
                            
                            if l==0:
                                if i ==0:
                                    R_l0n0[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l0n1[j][n]= Ri[j][n]
                                if i ==2:
                                    R_l0n2[j][n]= Ri[j][n]
                            if l==1:
                                if i ==0:
                                    R_l1n0[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l1n1[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l1n2[j][n]= Ri[j][n]
                            if l==2:
                                if i ==0:
                                    R_l2n0[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l2n1[j][n]= Ri[j][n]
                                if i ==1:
                                    R_l2n2[j][n]= Ri[j][n]
    
    Rl0 =[R_l0n0,R_l0n1,R_l0n2]
    Rl1 =[R_l1n0,R_l1n1,R_l1n2]
    Rl2 =[R_l2n0,R_l2n1,R_l2n2]
    R =[Rl0,Rl1,Rl2]
    print()
    for i in range(3):  
        print("Riemann selon",k[i],":")
        pprint(R[i])
        print()                 
    return R 

def Ricci(R):

    ricci1=[[0,0,0],
          [0,0,0],
          [0,0,0],]
    ricci=[[0,0,0],
          [0,0,0],
          [0,0,0],]

    """ Ricci """
    print("\u0332".join("Tenseur de Ricci :"))
    for l in range(3):
        for i in range(3):
            for j in range(3):
                ricci1[i][j] = R[l][i][j][l]
                if ricci1[i][j]!=0:
                    ricci[i][j]=ricci1[i][j]
                        
    print(ricci)
    print()
    return ricci

def Courbure_scalaire(ricci):
    scalaire=[[0,0,0],
              [0,0,0],
              [0,0,0],]
    
    """ courbure  scalaire """
    
    print("\u0332".join("Courbure scalaire :"))
    for i in range(3):
        for j in range(3):
            scalaire[i][j]= g_inv[i,j]* ricci[i][j]
            
    somme = scalaire[0][0]+scalaire[1][0]+scalaire[2][0]+scalaire[0][1]+scalaire[1][1]+scalaire[2][1]+scalaire[0][2]+scalaire[1][2]+scalaire[2][2]
    print(somme)
    return somme


"""--------------- Partie pour générer le code LaTeX ----------------- """

def PythonToLaTeX(Christofell,R,ricci,somme):
        
    def pmatrix(a):
        
        lines = str(a).replace('[', '').replace(']', '').replace(' - ','-').replace(' + ','+').splitlines()
        rv = [r' $$\begin{pmatrix}']
        rv += ['  ' + ' ,& '.join(l.split()) + r'\\' for l in lines]
        rv +=  [r'\end{pmatrix}$$']   
        #rv +=  [r'\end{pmatrix}$$'] 
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
    print('\section{RG3x3}',file=f)
    print('\subsection{Métrique}',file=f)
    print(pmatrix(np.array(g)),'\n',file=f)
    print('\subsection{Christofell}',file=f)
    for m in range(3):
        for i in range(3):
            for j in range(3):
                if Christofell[m][i][j] != 0:
                    print('$\Gamma^{',k[m],'}_{',k[i],k[j],'}',"=",Christofell[m][i][j],r'$\\', file=f)
                else :pass
        print('\n',file=f)
            
            
    #print('Symboles de Christofell',Tensor_matrix(np.array(Christofell)),file=f)
    for i in range(3):
        print('\paragraph{Symbole de Christofell selon $',k[i],'$ :}\n',pmatrix(np.array(Christofell[i])) + '\n', file=f)
    
    print('\subsection{Riemann}',file=f)
    for l in range(3):
        for i in range(3):
            for j in range(3):
                for n in range(3):
                    if R[l][i][j][n] !=0:
                        print('$R^{',k[l],'}_{',k[i],k[j],k[n],'}',"=",R[l][i][j][n],r'$\\', file=f)
        print('\n',file=f)
        
    #print('Tenseur de Riemann',Tensor_matrix(np.array(R)),file=f)
        
    for i in range(3):
        print('\paragraph{Tenseur de Riemann selon $',k[i],'$ :}\n',Tensor_matrix(np.array(R[i])) + '\n', file=f)
        
    print('\subsection{Ricci}',file=f)
    for l in range(3):
        for i in range(3):
            for j in range(3):
                if R[l][i][j][l] !=0:
                        print('$R\mathrm{icci}_{',k[i],k[j],'}',"=",R[l][i][j][l],r'$\\', file=f)
        print('\n',file=f)
    
    print('\paragraph{Tenseur de Ricci} ','\n',pmatrix(np.array(ricci)) + '\n', file=f)
    print('\paragraph{Courbure scalaire} ','\n','$$',somme , '$$',file=f)
    
        
    f.close()



Christofell=Christo3x3()
R = Riemann(Christofell)
ricci=  Ricci(R)
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

  