#----------------------------------------------------------------------------------------------------------------#    
#                                 IMPORT DELLE LIBRERIE UTILI
#-----------------------------------------------------------------------------------------------------------------

#  numpy > libreria numerica per la manipolazione di array e matrici e per il metodo solve_banded() usato per 
#          risolvere la equazione algebrica finale
#  matplotlib.pyplot > libreria grafica
#  time > per la stima del tempo di esecuzione dei vari metodi
#  sys > per facilitare l'input da terminale
#  sympy > per la visualizzazione testuale della funzione su terminale
 
import numpy as np 
import matplotlib.pyplot as plt
import time
import sys
import sympy as sy

#----------------------------------------------------------------------------------------------------------------#    
#                                 DEFINIZIONE DELLA CLASSE myFEM
#-----------------------------------------------------------------------------------------------------------------

class myFEM:

#
# ----------------------------------------------------------------------------------------- Attributi della classe   
#

#parametri di default
    ro = 7.88
    c = 0.437
    k = 0.836 
#    

#array utili
    risultato = []        
    mb = []
    kb = []
    M = np.empty(shape=(99,99))
    K = np.empty(shape=(99,99))
    xmax = 100
#

#
#--------------------------------------------------------------------------------------------- Metodi della classe
#

#----------------------------------------------------------------------------------------------------------------#    
#                                 COSTRUTTORE DELLA CLASSE myFEM
#-----------------------------------------------------------------------------------------------------------------

#COSTRUTTORE
    def __init__(self, nodi=11, xmin=0, xmax=100):
        
        self.nodi = nodi
        self.xmin = xmin
        self.xmax = xmax
        xvals = np.linspace(xmin, xmax, nodi)
        self.xvals = xvals
        self.h = xvals[1] - xvals[0]
        self.print_mesh()
        """
        print("\nDi default, creo una griglia spaziale (mesh) all'interno del dominio di definizione"
              " che ha lunghezza ",xmax," :\n\n ",xvals)
        print("\nCaratterizzata da : ")
        print("\t\tpunto iniziale : ",xmin,"\n\t\tpunto finale : ",xmax,"\n\t\tnumero di nodi : ",nodi,
              "\n\t\tampiezza di ciascun intervallo : ",self.h,"\n")       
        """      
#    

#----------------------------------------------------------------------------------------------------------------#    
#                                 INSERIMENTO PARAMETRI DEL PROBLEMA
#-----------------------------------------------------------------------------------------------------------------

#RIASSEGNAZIONE PARAMETRI
    def query_yes_no(self,question,vecchio_input,default="yes"):
    
        valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}

        if default is None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)

        while True:
            sys.stdout.write(question + prompt)
            choice = input().lower()
            if default is not None and choice == "":
                if (valid[default]):
                    nuovo_input = input("Inserisci nuovo valore : ")
                    return nuovo_input
                else:
                    return vecchio_input
            elif choice in valid:
                if (valid[choice]):
                    nuovo_input = input("Inserisci nuovo valore : ")
                    return nuovo_input
                else:
                    return vecchio_input
            else:
                sys.stdout.write("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")
#

#INSERIMENTO PARAMETRI
    def set_xvals(self):
        xvals = np.linspace(int(self.xmin), int(self.xmax), int(self.nodi))
        self.xvals = xvals
        h = xvals[1] - xvals[0]
        self.h = h

    def set_nodi(self,nodi):
        self.nodi = int(nodi)

    def set_xmin(self,xmin):
        self.xmin = int(xmin)

    def set_xmax(self,xmax):
        self.xmax = int(xmax)

    def set_ro(self,ro):
        self.ro = ro
        
    def set_c(self,c):
        self.c = c
     
    def set_k(self,k):
        self.k = k       
#

#FORZANTE    
    def rhs(self,xx,t):
        return 10**(-8)*t*xx*(100-xx)**2
#     

#----------------------------------------------------------------------------------------------------------------#    
#                                 CREAZIONE DELLA BASE DI POLINOMI
#-----------------------------------------------------------------------------------------------------------------

#BASE DI POLINOMI LINEARI A PEZZI    
    def myPiecewise_linear(self,xx,i):
        h = self.h
        p1 = (1/h)*(xx-(i-1)*h)
        p2 = -(1/h)*(xx-(i+1)*h)
        return np.piecewise(xx,
                            [xx<(i-1)*h, xx>(i+1)*h, ((i-1)*h<=xx) & (xx<(i*h)), (i*h<=xx) & (xx<=(i+1)*h)],
                            [0, 0, lambda x:p1, lambda x:p2])
#

#----------------------------------------------------------------------------------------------------------------#    
#                                 PRODOTTO BASE X RHS
#-----------------------------------------------------------------------------------------------------------------
        
#PRODOTTO BASE DI POLINOMI A PEZZI X RHS
    def product (self,xx,t,i):
        rhs = self.rhs
        myPiecewise_linear = self.myPiecewise_linear
        return rhs(xx,t)*myPiecewise_linear(xx,i)
#

#----------------------------------------------------------------------------------------------------------------#    
#                              STAMPA GRAFICI - ARRAY - MATRICI
#-----------------------------------------------------------------------------------------------------------------

#STAMPA GRIGLIA SPAZIALE
    def print_mesh(self):

        print("Griglia spaziale : \n\n ",self.xvals)
        print("\nCaratterizzata da : ")
        print("\t\tpunto iniziale : ",self.xmin,"\n\t\tpunto finale : ",self.xmax,"\n\t\tnumero di nodi : ",self.nodi,
              "\n\t\tampiezza di ciascun intervallo : ",self.h,"\n")  
#

#STAMPA SU TEMRINALE DELLA FORZANTE    
    def print_rhs(self):
        print("Valori della forzante per ognuno dei ",self.nodi," punti della griglia spaziale, in funzione di t :\n")
        xx = self.xvals
        ff = myFEM.rhs(self,xx,t=sy.symbols('t'))
        print(ff) 
        print("\nLa funzione corrispondente è la seguente. Se la si desidera cambiare, prego farlo manualmente" 
              " all'interno del metodo print_rhs di myFEM in quanto non so tradurre  l'input da simbolico a numerico\n ")
        print(myFEM.rhs(self,xx=sy.symbols('x'),t=sy.symbols('t')))
#    
  

#STAMPA BASE SU GRAFICO IN FORMATO DINAMICO 
    def stampaBase (self):

        start_time = time.process_time()

        h = self.h
        x = self.xvals
        print('x = ',x)
        print('h = ',h)
        xmin = self.xmin
        xmax = self.xmax
        nodi =self.nodi
        myPiecewise_linear = self.myPiecewise_linear

        vec = np.vectorize(myPiecewise_linear)
        extr = int(xmax/(nodi-1))
        
        #interaction on 
        plt.ion()
        
        for k in range(xmin+1,nodi-1):
            xvec = vec(x,k)
            print("creazione elemento numero : ",k)
            plt.plot(x,xvec)
            plt.title('Base di polinomi lineari a pezzi')
            plt.grid('True')
            plt.pause(0.1)    
  
        plt.pause(5)  
        plt.close('all')   

        elapsed_time = time.process_time() - start_time
        print("Tempo di esecuzione : ",elapsed_time," s")
#  

#STAMPA SU GRAFICO DINAMICO DEL PRODOTTO
    def print_product(self,t):
        
        start_time = time.process_time()

        xmin = self.xmin
        nodi = self.nodi
        x = self.xvals
        product = self.product
        pvec = np.vectorize(product) 

        for tempo in range(t):
            for i in range(xmin+1,nodi-1):
                vect = pvec(x,i,tempo)
                    #print("temp vale ",temp)
                    #print("vect vale ",vect)
                plt.plot(x,vect)
                plt.pause(0.1)
                plt.grid(True)
                plt.title("Prodotto forzante x base di polinomi")
        plt.pause(0.4)
        plt.close('all')

        elapsed_time = time.process_time() - start_time
        print("Tempo di esecuzione : ",elapsed_time," s")
# 

#----------------------------------------------------------------------------------------------------------------#    
#                                 COSTRUZIONE MATRICI TRIDIAGONALI
#-----------------------------------------------------------------------------------------------------------------

#DEFINIZIONE MATRICE TRIDIAGONALE  

    def tridiagMatrix(self,T,x,y,z,k1=-1, k2=0, k3=1):
        a = [x]*(T-abs(k1)); b = [y]*(T-abs(k2)); c = [z]*(T-abs(k3))
        return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)
#

#DEFINIZIONE MATRICE K DI RIGIDITA'

    def Kmatrix (self):
        h = self.h
        k = self.k    
        xmax = self.xmax
        nodi = self.nodi
        tridiagMatrix = self.tridiagMatrix

        K = tridiagMatrix(nodi-1,-(k/h),(2*k)/h,-(k/h)) 
        self.K = K 
        print("\nLa matrice K di rigidità, di forma ",np.shape(K)," è la seguente :\n ")                         
        print(K) 

        ak=np.empty(len(K))
        uk=np.empty(len(K))
        lk=np.empty(len(K))

        for i in range (len(K)):
            for j in range (len(K[i])):
                if (i==j):
                    ak[j]=K[i][j]
                    #print("M[i][j]=",M[i][j])
                if (i==j-1):
                    uk[j]=K[i][j]
                    #print("M[i][j]=",M[i][j])
                if (i==j+1):
                    lk[j]=K[i][j]
                    #print("M[i][j]=",M[i][j])

        kb = np.vstack((uk,ak,lk))
        self.kb = kb
        #print("kb=",kb,type(kb),np.shape(kb))
#

#DEFINIZIONE MATRICE M DI MASSA

    def Mmatrix(self):   
        h = self.h
        ro = self.ro
        c = self.c
        xmax = self.xmax
        nodi = self.nodi
        tridiagMatrix = self.tridiagMatrix 

        M=tridiagMatrix(nodi-1,(h*ro*c)/6,(2*h*ro*c)/3,(h*ro*c)/6)    
        self.M = M  
        print("\nLa matrice M di massa, di forma ",np.shape(M)," è la seguente :\n ")                         
        print(M)
        
        am = np.empty(len(M))
        um = np.empty(len(M))
        lm = np.empty(len(M))

        for i in range (len(M)):
            for j in range (len(M[i])):
                if (i==j):
                    am[j]=M[i][j]
                    #print("M[i][j]=",M[i][j])
                if (i==j-1):
                    um[j]=M[i][j]
                    #print("M[i][j]=",M[i][j])
                if (i==j+1):
                    lm[j]=M[i][j]
                    #print("M[i][j]=",M[i][j])

        mb = np.vstack((um,am,lm))
        self.mb = mb
        #print("mb=",mb,type(mb),np.shape(mb))    
#

#----------------------------------------------------------------------------------------------------------------#    
#                                 INTERPOLANTE FATTA DI POLINOMI LINEARI A PEZZI
#-----------------------------------------------------------------------------------------------------------------

#DEFINIZIONE DEL METODO RISOLUTIVO    
    def solve(self,t):
        start_time = time.process_time()

    #RICHIAMO LE FUNZIONI CHE DEFINISCONO LE MATRICI M K 
        self.Mmatrix()
        self.Kmatrix()

        K = self.K #3diag
        ab = self.mb 
    #

    #INVOCO GLI ATTRIBUTI GLOBALI DI CLASSE 
        nodi = self.nodi
        xmin = self.xmin
        xmax = self.xmax
        incr = int(xmax/(nodi-1))    
        x = self.xvals
        product = self.product
    # 

    #DEFINISCO GLI ARRAY UTILI        
        alpha = np.zeros(nodi-1)
        b = np.zeros(nodi-1)
        sarray = np.zeros(nodi-1)        
    #

    #INFO SUGLI STEP TEMPORALI
        dt = 0.69
        nt_step = int(t/dt)
        print("\n\nintervallo temporale minimo considerato = ",dt)
        print('tempo totale considerato = ',t)
        print("numeri di step nel tempo = ",nt_step) 
    #

    #RISOLVO LA EQUAZIONE ALGEBRICA ASSOCIATA AL PROBLEMA DIFFERENZIALE, CALCOLANDO I VALORI ASSUNTI DAI
    #COEFFICIENTI ALPHA        
        from scipy.integrate import quad
        from scipy.linalg import solve_banded
       
        f = np.zeros(nodi-1)

        for tempo in range(0,nt_step):
            for k in range(xmin+1,nodi-1): 
                index = int(k*incr)
                inf = index-incr
                sup = index+incr
                #print('inf = ',inf)
                #print('sup = ',sup)
                f[k-1] = quad(product,inf,sup,args=(tempo*dt,k))[0]
                #print('f dentro = ',f[k-1])
            print('f = ',f)

            b = -np.dot(K,alpha) + f
            sarray = solve_banded((1,1),ab,b)   
            alpha = alpha + dt*sarray           
    #
        
    #SOLUZIONE NUMERICA (INTERPOLANTE DI POLINOMI LINEARI)        
        def solution(self,xx):
            xmax = self.xmax
            nodi = self.nodi
            myPiecewise_linear = self.myPiecewise_linear
            sol = 0
            for i in range(0,nodi-1):
                sol = sol + alpha[i]*myPiecewise_linear(xx,i+1)
            return sol            
    #

    #STIMA DEL TEMPO DI ESECUZIONE
        elapsed_time = time.process_time() - start_time
        print("Tempo di esecuzione : ",elapsed_time," s")
    #    
#

#----------------------------------------------------------------------------------------------------------------#    
#                                 VISUALIZZAZIONE DELLA SOLUZIONE
#-----------------------------------------------------------------------------------------------------------------

    #STAMPA GRAFICO FINALE
        vecs = np.vectorize(solution)
        sol_vec = vecs(self,x)
        
        plt.plot(x,sol_vec)
        plt.xlabel('dominio spaziale')
        plt.ylabel('soluzione numerica')
        plt.title('Funzione di polinomi lineari interpolante')
        plt.grid(True)
        plt.pause(120)
    #      
#FINE

