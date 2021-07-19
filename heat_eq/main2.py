#-----------------------------------------------------------------------------------------------------------------#
                                                                                                                                 
# <  Risoluzione della equazione di propagazione del calore in 1 dimensione tramite FEM                                        
#   < Esercizio svolto nr. 6.8 della seconda edizione di Partial Differential Equations di Mark Gockenback                   
#      < Giada Costantini                                                                         
                                                                                                                  
#-----------------------------------------------------------------------------------------------------------------#

print("\nBenvenuto nel programma di risoluzione per problemi di propagazione del calore di un corpo rigido in"
      " una dimensione!\n\nIl programma è in grado di risolvere - applicando il metodo degli elementi finiti -"
      " problemi differenziali parziali al secondo ordine \nnon omogenei, espressi in FORMA FORTE e BEN POSTI, del tipo : "
      " \n\n\t\trc*du(x,t)/dt-k*d^2u(x,t)/dt^2 = f(x,t)\n")
      
print(" Con \nr = densità di materia [kg/m^3]\nc = calore specifico [J/K*kg] \nk = conducibilità termica [W/m*K]\n"
      " e \nu(x,t) = profilo di temperatura \nf(x,t) = sorgente di calore\n")

#----------------------------- struttura della equazione tipo ---------------------------------------------------#   
#   
#   rc* du(x,t)/dt - k* d^2u(x,t)/dt^2 = f(x,t)
#
#   I parametri del problema sono personalizzabili dall'utente. Di default prende i parametri dell'esercizio svolto
#   6.8 della seconda edizione del Gockenbach.   
#
# -----------------------------------------------------------------------------------------------------------------    

print("Il programma è composto da un file eseguibile indicato con main e un file di sola intestazione all'interno "
      "del quale \nvengono importate le librerie, definita la classe myFEM e i metodi per la risoluzione del problema.\n\n")

#---------------------------- struttura della cartella di lavoro -------------------------------------------------#    
#
#   (cartella)esame_python
#                   |_______> main2.py //file eseuguibile
#                   |_______> (cartella)moduli
#                                       |______> __init__.py //file vuoto
#                                       |______> myFEM2.py //header file
#                                       |______> (cartella)__pycache__ //creata di default
# 
# -----------------------------------------------------------------------------------------------------------------   


#----------------------------------------------------------------------------------------------------------------#    
#                                 IMPORTO IL FILE HEADER
#-----------------------------------------------------------------------------------------------------------------

from moduli import myFEM2

#----------------------------------------------------------------------------------------------------------------#    
#                                 ISTANZIO L'OGGETTO fem
#-----------------------------------------------------------------------------------------------------------------

fem = myFEM2.myFEM()
comments =  """
                  Di default l'oggetto fem costruisce un mesh spaziale (in 1 dimensione) composto da 11 nodi
                  che partono da 0 fino al valore 100, suddividendo il dominio in sotto-intervalli costanti 
                  di dimensione 10
            """

#----------------------------------------------------------------------------------------------------------------#    
#                        PERSONALIZZAZIONE E VISUALIZZAZIONE GRIGLIA SPAZIALE
#-----------------------------------------------------------------------------------------------------------------

#Visualizzo la griglia di default --------------------------------------------------------------------------------
fem.print_mesh()

#Ridefinisco, se desidero, i parametri che caratterizzano la griglia : nr nodi, xmin, xmax ----------------------
fem.set_nodi(fem.query_yes_no("Si desidera cambiare il numero dei nodi?",fem.nodi))
print("Va bene! Il nuovo valore dei nodi è : ",fem.nodi)

fem.set_xmin(fem.query_yes_no("Vuoi modificare xmin?",fem.xmin))
print("Va bene! Il nuovo valore di xmin è : ",fem.xmin)

fem.set_xmax(fem.query_yes_no("Vuoi modificare xmax?",fem.xmax))
print("Va bene! Il nuovo valore di xmax è : ",fem.xmax)

#Aggiorno il mesh spaziale, invocando il metodo di creazione della griglia --------------------------------------- 
fem.set_xvals()

#Ristampo su temrinale la nuova griglia --------------------------------------------------------------------------
fem.print_mesh()

#----------------------------------------------------------------------------------------------------------------#    
#                               PERSONALIZZAZIONE DEI PARAMETRI DEL PROBLEMA
#-----------------------------------------------------------------------------------------------------------------

#Aggiorno, se desidero, il valore della densità di materia -------------------------------------------------------
print("Il valore di default di ro è : ",fem.ro)
fem.set_ro(fem.query_yes_no("Vuoi aggiornare ro?",fem.ro))
print("ok. Il valore di ro è : ",fem.ro)

#Aggiorno, se desidero, il valore del calore specifico -----------------------------------------------------------
print("Il valore di default di c è : ",fem.c)
fem.set_c(fem.query_yes_no("Vuoi aggiornare c?",fem.c))
print("ok. Il valore di c è : ",fem.c)

#Aggiorno, se desidero, il valore della conducibilità termica ----------------------------------------------------
print("Il valore di default di k è : ",fem.k)
fem.set_k(fem.query_yes_no("Vuoi aggiornare k?",fem.k))
print("ok. Il valore di k è : ",fem.k)

#Visualizzazione della forzante
fem.print_rhs()

#----------------------------------------------------------------------------------------------------------------#    
#                   VISUALIZZAZIONE GRAFICA DINAMICA DELLA BASE DI POLINOMI LINEARI A PEZZI
#-----------------------------------------------------------------------------------------------------------------

fem.stampaBase()  
comments = """
                Visualizzo la base di polinomi lineari a pezzi. La base è stata creata usando il metodo 
                piecewise della libreria numpy e le funzioni lambda. La base non dipende dal parametro temporale.
           """

#----------------------------------------------------------------------------------------------------------------#    
#                     CREAZIONE E STAMPA SU TERMINALE DELLE MATRICI DI MASSA E RIGIDITà
#-----------------------------------------------------------------------------------------------------------------

fem.Kmatrix()
fem.Mmatrix()
comments = """
                Le matrici K ed M sono tridiagonali poiché si è scelta una base triangolare del dominio.
           """

#----------------------------------------------------------------------------------------------------------------#    
#                              RISOLUZIONE NUMERICA E STAMPA GRAFICO 
#-----------------------------------------------------------------------------------------------------------------

fem.solve(180)
comments = """
                Invoco il metodo risolutivo che richiede un solo argomento: il tempo a cui si vuole valutare 
                la distribuzione del calore nel corpo. 
           """


