MODULE save_file
!Subroutine per la scrittura su file
IMPLICIT NONE

    CONTAINS

    SUBROUTINE export_mat(mat, T_mat, dim1, dim2, filename, opt) !Si occupa della scrittura dei vari file di output
        INTEGER, INTENT(IN) :: dim1, dim2, opt
        REAL*8, INTENT(IN) :: mat(dim1,dim2), T_mat(dim1)
        CHARACTER(LEN=*), INTENT(IN) :: filename

        INTEGER :: i, j
        CHARACTER(LEN=4), PARAMETER :: extension = '.dat'
        CHARACTER(LEN=20), DIMENSION(dim2+1) :: labels

        SELECT CASE(opt)

            CASE(0)
                labels = ['Log10(T)            ', 'nH0/nH_tot          ', 'nH+/nH_tot          ', 'nHe0/nHe_tot        ',&
                        'nHe+/nHe_tot        ', 'nHe++/nHe_tot       ']

            CASE(1)
                labels = ['Log10(T)            ', 'Log10(coll_exc_H0)  ', 'Log10(coll_exc_He+) ', 'Log10(coll_ion_H0)  ',&
                        'Log10(coll_ion_He0) ', 'Log10(coll_ion_He+) ', 'Log10(recomb_H+)    ', 'Log10(recomb_He+)   ',&
                        'Log10(recomb_He++)  ', 'Log10(diel_rec_He+) ', 'Log10(freefree_ions)', 'Log10(cooling_tot)  ']

            CASE(2)
                labels = ['Log10(T)            ', 'Log10(coll_exc_H0)  ', 'Log10(coll_ion_H0)  ', 'Log10(recomb_H+)    ',&
                        'Log10(freefree_ions)', 'Log10(cooling_tot)  ']

            CASE(3)
                labels = ['Log10(T)            ', 'Log10(cooling_time) ']

        END SELECT

        OPEN(11, file=filename//extension)
        WRITE(11,'(1x, *(g0, ", "))') (labels(i), i=1,SIZE(labels))
        DO i = 1, SIZE(mat,1)
            WRITE(11,'(1x, *(g0, ", "))') T_mat(i), (mat(i,j), j=1,SIZE(mat,2))    
        END DO
        CLOSE(11)


    END SUBROUTINE export_mat
    
END MODULE save_file

MODULE sistem_res
!Insieme di subroutine necessarie per calcolare lo stato di ionizzazione e la curva di cooling
IMPLICIT NONE

    CONTAINS

    SUBROUTINE check_mat(a,c,dim,order) !Controlla e riordina la matrice 'a' ed il vettore 'c' affinché siano compatibili con l'eliminazione di Gauss
        INTEGER, INTENT(IN) :: dim
        REAL*8, INTENT(INOUT) :: a(dim,dim), c(dim)
        INTEGER, INTENT(OUT) :: order(dim) !'order' è un vettore che, una volta esportato, permette di risalire all'ordine originale del vettore delle soluzioni

        INTEGER :: i, j, pos, k
        INTEGER :: temp1(dim), min_val
        REAL*8, DIMENSION(dim+1) :: temp2

        order = 0

        temp1 = 0
        DO i = 1, dim !Controlla la posizione del primo valore non nullo per ogni riga della matrice
            DO j = 1, dim
                IF (a(i,j) /= 0) THEN 
                    temp1(i) = j
                    EXIT
                END IF
            END DO
            order(i) = i
        END DO

        DO i = 1, dim-1 !Riordina la matrice 'a' ed il vettore 'c' seguendo il vettore 'temp1' in ordine crescente
            min_val = temp1(i)
            temp2(1:dim) = a(i,:)
            temp2(dim+1) = c(i)
            k = order(i)
            pos = i
            DO j = i+1, dim
                IF (temp1(j) .lt. min_val) THEN
                    min_val = temp1(j)
                    temp2(1:dim) = a(j,:)
                    temp2(dim+1) = c(j)
                    k = order(j)
                    pos = j
                END IF
            END DO

            temp1(pos) = temp1(i)
            a(pos,:) = a(i,:)
            c(pos) = c(i)
            order(pos) = order(i)

            temp1(i) = min_val
            a(i,:) = temp2(1:dim)
            c(i) = temp2(dim+1)
            order(i) = k

        END DO

    END SUBROUTINE check_mat

    SUBROUTINE gauss(a,c,x,ndati) !Risolve con Gauss il sistema richiesto
        INTEGER :: ndati, i, j, k
        REAL*8, DIMENSION(ndati) :: c, x !'c' è il vettore dei termini noti e 'x' è il vettore delle soluzioni
        REAL*8, DIMENSION(ndati,ndati) :: a !'a' è la matrice del sistema
        REAL*8 :: fakt, summa

        DO i = 1,ndati-1 !Sceglie la variabile da eliminare
            DO j = i+1,ndati !Sceglie la riga su cui eliminare
                fakt = a(j,i)/a(i,i)
                DO k = 1,ndati
                    a(j,k) = a(j,k) - a(i,k)*fakt
                END DO
                c(j) = c(j) - c(i)*fakt
            END DO
        END DO
    
        x(ndati) = c(ndati)/a(ndati,ndati)
        DO i = ndati-1,1,-1
            summa = 0.d0
            DO j = i+1,ndati
                summa = summa + a(i,j)*x(j)
            END DO
            x(i) = (c(i) - summa)/a(i,i)
        END DO
    
    END SUBROUTINE gauss

    SUBROUTINE recombination_rates(T, rates_mat) !In funzione della temperatura calcola i tassi di ricombinazione
        REAL*8, INTENT(IN) :: T
        REAL*8, INTENT(OUT) :: rates_mat(7)

        REAL*8 :: par1, par2
        REAL*8 :: alpha_H1, alpha_He1, alpha_d, alpha_He2, Gamma_eH0, Gamma_eHe0, Gamma_eHe1 !Tassi di ricombinazione
    
        par1 = T**(-0.5)*(T*1.d-3)**(-0.2)*1.d0/(1.d0 + (T*1.d-6)**(0.7))
        par2 = 1.d0/(1.d0 + (T*1.d-5)**(0.5))
        
        alpha_H1 = 8.4d-11*par1
        alpha_He1 = 1.5d-10*T**(-0.6353)
        alpha_d = 1.9d-3*T**(-1.5)*EXP(-470000.0/T)*(1.d0 + 0.3d0*EXP(-94000.0/T))
        alpha_He2 = 3.36d-10*par1
        Gamma_eH0 = 5.85d-11*T**(0.5)*EXP(-157809.1/T)*par2
        Gamma_eHe0 = 2.38d-11*T**(0.5)*EXP(-285335.4/T)*par2
        Gamma_eHe1 = 5.68d-12*T**(0.5)*EXP(-631515.0/T)*par2
    
        rates_mat = [alpha_H1, alpha_He1, alpha_d, alpha_He2, Gamma_eH0, Gamma_eHe0, Gamma_eHe1]
    END SUBROUTINE recombination_rates

    SUBROUTINE fraction(rec_rates, dim1, res, dim2) !Calcola le frazioni degli stati di ionizzazione
        INTEGER, INTENT(IN) :: dim1, dim2
        REAL*8, DIMENSION(dim1), INTENT(IN) :: rec_rates
        REAL*8, DIMENSION(dim2), INTENT(OUT) :: res

        INTEGER :: i

        DO i = 1, SIZE(res)
            IF (i+1 .le. 3) THEN
                res(i) = rec_rates(i+1)/(rec_rates(2) + rec_rates(3))
            ELSE IF (i+1 .gt. 3 .and. i+1 .le. 6) THEN
                res(i) = rec_rates(i+1)/(rec_rates(4) + rec_rates(5) + rec_rates(6))    
            END IF
        END DO

    END SUBROUTINE fraction

    SUBROUTINE cooling_rates(T, n_mat, dim1, rates_mat) !In funzione della temperatura e degli stati di ionizzazione calcola i tassi di cooling
        INTEGER, INTENT(IN) :: dim1
        REAL*8, DIMENSION(dim1), INTENT(IN) :: n_mat
        REAL*8, INTENT(IN) :: T
        REAL*8, dimension(11), INTENT(OUT) :: rates_mat

        REAL*8 :: temp(10)
        REAL*8 :: par2, par1, gff
        REAL*8 :: cexc_H0, cexc_He1, cion_H0, cion_He0, cion_He1, rec_H1, rec_He1, rec_He2, drec_He1, ff_ions !Tassi di cooling

        par1 = T**(0.5)*(T*1.d-3)**(-0.2)*1.d0/(1.d0 + (T*1.d-6)**(0.7))
        par2 = 1.d0/(1.d0 + (T*1.d-5)**(0.5))
        gff = 1.1d0 + 0.34d0*EXP(-((5.5d0 - LOG10(T))**2)/3.d0)

        cexc_H0 = 7.50d-19*EXP(-118348.0/T)*par2*n_mat(1)*n_mat(2)
        cexc_He1 = 5.54d-17*T**(-0.397)*EXP(-473638.0/T)*par2*n_mat(1)*n_mat(5)
        cion_H0 = 1.27d-21*T**(0.5)*EXP(-157809.1/T)*par2*n_mat(1)*n_mat(2)
        cion_He0 = 9.38d-22*T**(0.5)*EXP(-285335.4/T)*par2*n_mat(1)*n_mat(4)
        cion_He1 = 4.95d-22*T**(0.5)*EXP(-631515.0/T)*par2*n_mat(1)*n_mat(5)
        rec_H1 = 8.70d-27*par1*n_mat(1)*n_mat(3)
        rec_He1 = 1.55d-26*T**(0.3647)*n_mat(1)*n_mat(5)
        rec_He2 = 3.48d-26*par1*n_mat(1)*n_mat(6)
        drec_He1 = 1.24d-13*T**(-1.5)*EXP(-470000.0/T)*(1.d0 + 0.3d0*EXP(-94000.0/T))*n_mat(1)*n_mat(5)
        ff_ions = 1.42d-27*gff*T**(0.5)*(n_mat(3) + n_mat(5) + 4.d0*n_mat(6))*n_mat(1)

        temp = [cexc_H0, cexc_He1, cion_H0, cion_He0, cion_He1, rec_H1, rec_He1, rec_He2, drec_He1, ff_ions]
        rates_mat = [temp, SUM(temp)]
        
    END SUBROUTINE cooling_rates

    SUBROUTINE buildmat(rec_rates_mat, size_mat, X, nH, n, size_n) !Costruisce la matrice necessaria per risolvere il sistema degli stati di ionizzazione
        INTEGER, INTENT(IN) :: size_mat, size_n
        REAL*8, INTENT(IN) :: X, nH
        REAL*8, DIMENSION(size_mat), INTENT(IN) :: rec_rates_mat
        REAL*8, DIMENSION(size_n), INTENT(OUT) :: n

        INTEGER :: order(size_n), i
        REAL*8 :: a(size_n,size_n), c(size_n), sol(size_n) !'a' è la matrice del sistema, 'c' è il vettore dei termini noti e 'sol' è il vettore delle soluzioni
        REAL*8 :: y

        y = (1.d0 - X)/(4.d0*X)
        a = 0.d0
        c = 0.d0
        c(4) = 1.d0
        c(6) = y

        !Definizione della matrice corrispondente al sistema richiesto
        a(1,:) = [0.d0, rec_rates_mat(5), -rec_rates_mat(1), 0.d0, 0.d0, 0.d0]
        a(2,:) = [0.d0, 0.d0, 0.d0, rec_rates_mat(6), -rec_rates_mat(2) - rec_rates_mat(3), 0.d0]
        a(3,:) = [0.d0, 0.d0, 0.d0, 0.d0, -rec_rates_mat(7), rec_rates_mat(4)]
        a(4,:) = [0.d0, c(4), c(4), 0.d0, 0.d0, 0.d0]
        a(5,:) = [c(4), 0.d0, -c(4), 0.d0, -c(4), -2*c(4)]
        a(6,:) = [0.d0, 0.d0, 0.d0, c(4), c(4), c(4)]

        CALL check_mat(a, c, size_n, order) !Riordina la matrice in modo che sia risolvibile con Gauss
        CALL gauss(a, c, sol, size_n) !Risolve con Gauss il sistema

        DO i = 1, SIZE(order) !Riordina il vettore delle soluzioni e lo moltiplica per nH
            n(order(i)) = nH*sol(order(i))
        END DO

        !n = [ne, nH0, nH1, nHe0, nHe1, nHe2]

    END SUBROUTINE buildmat

END MODULE sistem_res

MODULE ODE
!Insieme di funzioni e subroutine necessarie per risolvere l'equazione differenziale ordinaria
USE sistem_res
IMPLICIT NONE
    REAL*8, PRIVATE, PARAMETER :: gamma = 5.d0/3.d0 !Parametro gamma
    REAL*8, PRIVATE, PARAMETER :: mp = 1.67262192d-27 !Massa del protone
    REAL*8, PRIVATE, PARAMETER :: kB = 1.380649d-23 !Costante di Boltzmann

    CONTAINS

    REAL*8 FUNCTION mu_fun(X, nH, ne) !Calcola il peso molecolare medio
        REAL*8, INTENT(IN) :: X, nH, ne
        REAL*8 :: y

        y = (1.d0 - X)/(4.d0*X)
        mu_fun = (1.d0 + 4.d0*y)/(1.d0 + y + (ne/nH))

    END FUNCTION mu_fun

    REAL*8 FUNCTION u_fun(T, mu) !Calcola l'energia interna per unità di massa
        REAL*8, INTENT(IN) :: T, mu

        u_fun = 1.d0/(gamma - 1.d0)*(kB*T/(mu*mp))

    END FUNCTION u_fun

    REAL*8 FUNCTION rho_fun(nH, X) !Calcola la densità di massa
        REAL*8, INTENT(IN) :: nH, X

        rho_fun = nH*mp/X

    END FUNCTION rho_fun

    REAL*8 FUNCTION tcool_fun(u, rho, cooling) !Calcola il tempo di cooling
        REAL*8, INTENT(IN) :: u, rho, cooling

        tcool_fun = u*rho/cooling

    END FUNCTION tcool_fun

    REAL*8 FUNCTION fun(T, u, mu) !Funzione di cui bisogna trovare lo zero per ottenere la temperatura
        REAL*8, INTENT(IN) :: T, u, mu

        fun = T - (gamma - 1.d0)*(u*mu*mp)/kB

    END FUNCTION fun

    SUBROUTINE secant(T, u, X, nH) !Implementazione del metodo della secante per trovare lo zero di una funzione
        REAL*8, INTENT(IN) :: u, nH, X
        REAL*8, INTENT(INOUT) :: T

        INTEGER, PARAMETER :: max_iter = 1000 !Numero massimo di iterazioni consentite
        INTEGER :: iter
        REAL*8, PARAMETER :: h = 0.01d0, toll = 1.d0 !'h' è il passo con cui vene calcolato l'incremento e 'toll' è la tolleranza per determinare se il metodo è riuscito a convergere
        REAL*8 :: rates(7), n_arr(6)
        REAL*8 :: t1, t2, fun1, fun2, mu, var

        iter = 0
        var = toll + 1 !Si inizializza 'var' con un valore sicuramente diverso da 'toll'
        DO WHILE(var .gt. toll) !Il ciclo continua finché la variazione è minore o uguale alla tolleranza
            IF (iter == max_iter) THEN !Se si raggiungono il numero massimo di iterazioni il ciclo viene interrotto
                PRINT *, "Il programma non e' riuscito ad arrivare a convergenza."
                EXIT
            END IF

            !Si assume come temperatura quella iniziale e si calcola la funzione di cui bisogna ricercare lo zero
            t1 = T
            CALL recombination_rates(t1, rates)
            CALL buildmat(rates, SIZE(rates), X, nH, n_arr, SIZE(n_arr))
            mu = mu_fun(X, nH, n_arr(1))
            fun1 = fun(t1, u, mu)

            !Si aumenta la temperatura di un passo 'h' e si calcola nuovamente la funzione di cui bisogna ricercare lo zero
            t2 = t1 + h
            CALL recombination_rates(t2, rates)
            CALL buildmat(rates, SIZE(rates), X, nH, n_arr, SIZE(n_arr))
            mu = mu_fun(X, nH, n_arr(1))
            fun2 = fun(t2, u, mu)

            T = t1 - fun1*h/(fun2-fun1) !Si applica la formula moodificata del metodo della secante
            var = ABS(T-t1) !Si calcola la variazione tra la temperatura ottenuta con il metodo della secante e la temperatura iniziale
            
            iter = iter + 1
        END DO

    END SUBROUTINE secant

    SUBROUTINE rk2(T, time_val, T_min, u_min, t0, u0, cost, X, nH) !Implementazione del metodo di Runge-Kutta al secondo ordine
        REAL*8, INTENT(IN) :: T_min, u_min, t0, u0, cost, X, nH
        REAL*8, INTENT(INOUT) :: T, time_val

        REAL*8, PARAMETER :: b2 = 0.5d0 !Paramentro da cui dipendono gli altri parametri del metodo. Con 'b2 = 0.5d0' il metodo coincide con il metodo di Heun
        REAL*8, PARAMETER :: b1 = 1.d0 - b2, a11 = 1.d0/(2.d0*b2), c1 = 1.d0/(2.d0*b2)
        REAL*8 :: rates(7), n_arr(6), cooling(11)
        REAL*8 :: k1, k2, h, rho, mu, u_temp, t_cool, u

        CALL recombination_rates(T, rates)
        CALL buildmat(rates, SIZE(rates), X, nH, n_arr, SIZE(n_arr))
        CALL cooling_rates(T, n_arr, SIZE(n_arr), cooling)
        
        rho = rho_fun(nH, X)
        mu = mu_fun(X, nH, n_arr(1))
        u_temp = u_fun(T, mu)
        t_cool = tcool_fun(u_temp, rho, cooling(11))

        k1 = cost*cooling(11)/(nH**2)
        h = (t_cool/t0)*1.d-3 !Per la scelta del time step viene considerato il tempo di cooling adimensionalizzato

        time_val = time_val + c1*h
        u = u_temp/u0 + a11*h*k1 !Per il calcolo dell'energia interna viene considerata l'energia interna, calcolata inizialmente, adimensionalizzata
        CALL secant(T, u*u0, X, nH) !La temperatura viene aggiornata con il metodo della secante

        T = MAX(T, T_min)
        IF (T == T_min) THEN !Se la temperatura è pari alla temperatura minima allora l'energia interna è pari all'energia minima
            u = u_min
        ELSE
            CALL recombination_rates(T, rates)
            CALL buildmat(rates, SIZE(rates), X, nH, n_arr, SIZE(n_arr))
            CALL cooling_rates(T, n_arr, SIZE(n_arr), cooling)

            k2 = cost*cooling(11)/(nH**2)
            u = u_temp/u0 + h*(b1*k1 + b2*k2) !Con rk2 viene aggiornata l'energia interna, sempre considerando l'energia interna adimensionalizzata
            CALL secant(T, u*u0, X, nH) !La temperatura viene aggiornata con il metodo della secante

            T = MAX(T, T_min)
            IF (T == T_min) u = u_min !Se la temperatura è pari alla temperatura minima allora l'energia interna è pari all'energia minima
        END IF

    END SUBROUTINE rk2

    SUBROUTINE ODE_solver(params, arr, max_iter)
        REAL*8, INTENT(IN), DIMENSION(4) :: params
        INTEGER, INTENT(IN) :: max_iter
        REAL*8, INTENT(OUT) :: arr(max_iter,2)

        INTEGER :: iter
        REAL*8 :: rates(7), n_arr(6), cooling(11)
        REAL*8 :: u0, u_min, t0, mu, rho, cost, T, time_val
        REAL*8 :: T_min, T_max, X, nH !Valori che vengono passati come parametri della subroutine

        T_min = params(1)
        T_max = params(2)
        X = params(3)
        nH = params(4)
        arr = -99.9d0

        !Si calcola l'energia interna per la temperatura minima
        CALL recombination_rates(T_min, rates)
        CALL buildmat(rates, SIZE(rates), X, nH, n_arr, SIZE(n_arr))
        mu = mu_fun(X, nH, n_arr(1))
        u_min = u_fun(T_min, mu)

        !A temperatura massima si calcola l'energia interna, il tempo di cooling e la constante moltiplicativa per adimensionalizzare l'equazione differenziale
        CALL recombination_rates(T_max, rates)
        CALL buildmat(rates, SIZE(rates), X, nH, n_arr, SIZE(n_arr))
        CALL cooling_rates(T_max, n_arr, SIZE(n_arr), cooling)
        mu = mu_fun(X, nH, n_arr(1))
        u0 = u_fun(T_max, mu)
        rho = rho_fun(nH, X)
        t0 = tcool_fun(u0, rho, cooling(11))
        cost = -(t0/u0)*((nH**2)/rho)

        T = T_max
        iter = 0
        time_val = 0.d0
        DO WHILE (T .gt. T_min) !Il ciclo continua finché la temperatura non raggiunge la temperatura minima
            IF (iter == max_iter) THEN !Se si raggiunge il numero massimo di iterazioni il ciclo viene interrotto
                PRINT *, "Il programma ha raggiunto il numero massimo di iterazioni."
                EXIT
            END IF

            CALL rk2(T, time_val, T_min, u_min, t0, u0, cost, X, nH) !Ad ogni ciclo calcola la temperatura ed il tempo
            !Vengono salvati i risultati di rk2 in scala logaritmica
            arr(iter+1,1) = LOG10(T)
            arr(iter+1,2) = LOG10(time_val)

            iter = iter + 1
        END DO

    END SUBROUTINE ODE_solver

END MODULE ODE

PROGRAM main
USE save_file
USE sistem_res
USE ODE
IMPLICIT NONE
    INTEGER, PARAMETER :: npoints = 200 !Numero di punti in cui calcolare lo stato di ionizzazione e la curva di cooling
    INTEGER, PARAMETER :: max_iter = 7700 !Numero di iterazioni massime per il calcolo dell'evoluzione della temperatura del gas
    INTEGER :: i, j, k

    REAL*8, PARAMETER :: T_max = 1.d8, T_min = 1.d4 !Temperatura massima e minima per il calcolo dello stato di ionizzazione e della curva di cooling
    REAL*8, ALLOCATABLE :: temp(:,:), evolution_mat(:,:), ODE_params(:)
    REAL*8 :: rec_rates_mat076(0:npoints, 7), n076(0:npoints, 6), frac_mat076(0:npoints, 5)
    REAL*8 :: cool_rates_mat076(0:npoints, 11), T_mat(0:npoints)
    REAL*8 :: rec_rates_mat100(0:npoints, 7), n100(0:npoints, 6)
    REAL*8 :: cool_rates_mat100(0:npoints, 5)
    REAL*8 :: T, step, X, nH

    CHARACTER(LEN=28), DIMENSION(4) :: filenames

    ALLOCATE(temp(0:npoints,11))
    !Inizializza le varie matrici
    T_mat = -99.9d0
    n076 = -99.9d0
    frac_mat076 = -99.9d0
    rec_rates_mat076 = -99.9d0
    cool_rates_mat076 = -99.9d0

    n100 = -99.9d0
    rec_rates_mat100 = -99.9d0
    cool_rates_mat100 = -99.9d0

    step = (LOG10(T_max) - LOG10(T_min))/npoints !Suddivide l'intervallo tra la temperatura massima e la temperatura minima in modo da avere 'npoints' punti equidistanziati in scala logaritmica
    T = T_min
    nH = 1.d0

    PRINT *, "Inizio del calcolo dello stato di ionizzazione di 'H' e 'He' e delle curve di cooling in corso..."
    PRINT *, "Verranno eseguiti i calcoli per temperature comprese tra", T_min, "K e", T_max, "K"
    DO i = 0, npoints
        T = 10**(4+i*step) !Aggiorna la temperatura con cui fare i calcoli
        T_mat(i) = LOG10(T)

        X = 0.76d0 !Calcola lo stato di ionizzazione e le curve di cooling per il caso X = 0.76
        CALL recombination_rates(T, rec_rates_mat076(i,:))
        CALL buildmat(rec_rates_mat076(i,:), SIZE(rec_rates_mat076,2), X, nH, n076(i,:), SIZE(n076,2))
        CALL fraction(n076(i,:), SIZE(n076,2), frac_mat076(i,:), SIZE(frac_mat076,2))
        CALL cooling_rates(T, n076(i,:), SIZE(n076,2), cool_rates_mat076(i,:))

        X = 1.00d0 !Calcola le curve di cooling per il caso X = 1
        CALL recombination_rates(T, rec_rates_mat100(i,:))
        CALL buildmat(rec_rates_mat100(i,:), SIZE(rec_rates_mat100,2), X, nH, n100(i,:), SIZE(n100,2))
        CALL cooling_rates(T, n100(i,:), SIZE(n100,2), temp(i,:))

        k = 1
        DO j = 1, SIZE(temp,2)
            IF (temp(i,j) /= 0.d0) THEN !Per X = 1 salva solamente i tassi di cooling che non sono esattamente nulli
                cool_rates_mat100(i,k) = temp(i,j)
                k = k + 1
            END IF
        END DO

    END DO

    PRINT *, "Operazione completata con successo."
    PRINT *, ""
    !Salva su file i dati calcolati per lo stato di ionizzazione e le curve di cooling
    CALL export_mat(frac_mat076, T_mat, SIZE(frac_mat076,1), SIZE(frac_mat076,2), 'rapporti_nH010_X076', 0)
    PRINT *, "I dati dello stato di ionizzazione per X = 0.76 sono stati salvati nel file 'rapporti_nH010_X076.dat'"
    CALL export_mat(LOG10(cool_rates_mat076), T_mat, SIZE(cool_rates_mat076,1), SIZE(cool_rates_mat076,2), 'cooling_nH010_X076', 1)
    PRINT *, "I dati della curva di cooling per X = 0.76 sono stati salvati nel file 'cooling_nH010_X076.dat'"
    CALL export_mat(LOG10(cool_rates_mat100), T_mat, SIZE(cool_rates_mat100,1), SIZE(cool_rates_mat100,2), 'cooling_nH010_X100', 2)
    PRINT *, "I dati della curva di cooling per X = 1.00 sono stati salvati nel file 'cooling_nH010_X100.dat'"

    ALLOCATE(evolution_mat(max_iter,8), ODE_params(4))
    PRINT *, ""
    PRINT *, "Sono state impostate", max_iter, "iterazioni come numero massimo di iterazioni per risolvere le seguenti ODE."
    PRINT *, ""
    filenames = ['evoluzione_nH001_Tin1e6_X076', 'evoluzione_nH010_Tin1e6_X076', 'evoluzione_nH100_Tin1e6_X076', &
                'evoluzione_nH010_Tin1e7_X076'] !Array dei nomi con cui verranno salvati i dati calcolati dell'ODE

    PRINT *, "Inizio del calcolo dell'evoluzione della temperatura del gas in funzione del tempo in corso..."
    DO i = 1, 8, 2
        SELECT CASE (i) !Seleziona i parametri con cui calcolare l'ODE
            CASE(1)
                ODE_params = [1.d4, 1.d6, 0.76d0, 0.1d0]
            CASE(3)
                ODE_params = [1.d4, 1.d6, 0.76d0, 1.d0]
            CASE(5)
                ODE_params = [1.d4, 1.d6, 0.76d0, 10.0d0]
            CASE(7)
                ODE_params = [1.d4, 1.d7, 0.76d0, 1.d0]
        END SELECT

        PRINT *, "Dati della composizione chimica del gas:"
        PRINT *, "nH =", ODE_params(4), "X =", ODE_params(3)
        PRINT *, "e con temperature iniziale e finale rispettivamente:"
        PRINT *, "T_in =", ODE_params(2), "T_fin = ", ODE_params(1) 
        CALL ODE_solver(ODE_params, evolution_mat(:,i:i+1), max_iter)
        CALL export_mat(evolution_mat(:,i+1), evolution_mat(:,i), SIZE(evolution_mat,1), 1, filenames(NINT(i/2.d0)), 3) !Salva su file i calcoli appena effettuati
        PRINT *, "I dati dell'evoluzione temporale della temperatura sono stati salvati nel file:"
        PRINT *, "      '"//filenames(NINT(i/2.d0))//".dat'"
        PRINT *, ""
    END DO
    PRINT *, "Operazione completata con successo"

END PROGRAM main