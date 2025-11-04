#Codigo modelacion tanque de hidrogeno
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

#Parametros del tanque (preguntar bien si estos estan definidos o como)
FLUIDO = 'Hydrogen'
V =    111        #(m^3)
AREA_S =    111 #(m^2)
U_TW = 111  #(W/m^2*K)
T_AMB =     111  #(K)

#Parametros iniciales (preguntarle a la maite si estos son los que me da ella)
T0 = 111     #(K)
P0 = 111     #(Pa)  
y0 = [T0, P0] #Vector con las propiedades iniciales
T_IN= 111  #(K)
P_IN= 111   #(Pa)
H_IN= PropsSI('Hmass', 'T', T_IN, 'P', P_IN, FLUIDO) #(J/kg)
M_IN = 111  #(kg/s)  
M_OUT = 111 #(kg/s)


#Definimos el sistema de Ecuaciones Diferenciales Ordinarias (EDO) 
def modelo_tanque(t, y, V, U_tw, A_s, T_amb, h_in, m_in_cte, m_out_cte, fluid):
    
    
    #Explicacion de los parametros de entrada:
    #t: tiempo (s)
    #y: vector de las condiciones iniciales [T, P]
    #V: volumen del tanque (m^3)
    #U_tw: coef. de transferencia de calor (W/m^2*K)
    #A_s: area superficial del tanque (m^2)
    #T_amb: temperatura ambiente (K)
    #h_in: entalpía de entrada (J/kg)
    #m_in_cte: masa que esta entrando en el tiempo t
    #m_out_cte: masa que esta saliendo en el tiempo t
    #fluid: nombre del fluido para CoolProp
    
    
    #Definimos T y P desde el vector y
    T, P = y[0], y[1]
    
    #Aca deberia definir m_in y m_out, revisar como deben ser definidas (por mientras los dejo como una constante)
    m_in = m_in_cte
    m_out = m_out_cte
    
    #Vemos las propiedades termodinámicas en (T, P) con CoolProp
    try:
        #Propiedades
        rho = PropsSI('Dmass', 'T', T, 'P', P, fluid) #Densidad (kg/m^3)
        h = PropsSI('Hmass', 'T', T, 'P', P, fluid)   #Entalpia (J/kg)
        Cp = PropsSI('Cpmass', 'T', T, 'P', P, fluid) #Capacidad calorifica (J/kg*K)
        
        #Derivadas parciales
        drho_dT = PropsSI('d(Dmass)/d(T)|P', 'T', T, 'P', P, fluid)
        drho_dP = PropsSI('d(Dmass)/d(P)|T', 'T', T, 'P', P, fluid)
        dh_dP_T = PropsSI('d(Hmass)/d(P)|T', 'T', T, 'P', P, fluid)

    #Para definir casos donde CoolProp no funcione    
    except ValueError:
        print(f"Error de CoolProp en T={T:.2f} K, P={P:.2f} Pa. Deteniendo.")
        return [0, 0] #derivadas nuals

    #Para resolver el sistema se define A*x = b 
    #Construccion de la matriz de las ecuaciones diferenciales -> Izquierda de las ecuaciones (lo que multiplica dP/dt y dT/dt)
    
    #Balance de masa: V * ( drho/dT_P *dT/dt + drho/dP_T*dP/dt) = m_in - m_out
    A11 = V * drho_dT
    A12 = V * drho_dP
    
    #Balance de energia: rho*V*(Cp*dT/dt + dh/dP_T * dP/dt) - V*dP/dt = m_in*(h_in - h) - U_tw*A_s*(T - T_amb)
    A21 = rho * V * Cp
    A22 = (rho * V * dh_dP_T) - V
    
    #Matriz A
    A = np.array([
        [A11, A12],
        [A21, A22]
    ])
    
    
    #Ahora el vector b (lado derecho de las ecuaciones)
    #Blance de masa: 
    b1 = m_in - m_out
    
    #Balance de nergia:
    Q_loss = U_tw * A_s * (T - T_amb)
    
    b2 = m_in * (h_in - h) - Q_loss
    
    #Vector b
    b = np.array([b1, b2])
    
    #Ahora para resolver el sistema con respecto a dT/dt y dP/dt
    try:
        [dTdt, dPdt] = np.linalg.solve(A, b)
        return [dTdt, dPdt]
    #En caso de error pq la matriz es singular 
    except np.linalg.LinAlgError:
        print(f"Error: Matriz singular en T={T:.2f} K, P={P:.2f} Pa.")
        return [0, 0]
    


#En caso de que el flujo de entrada y salida sean funciones del tiempo o algo asi
def flujo_entrada_func(t):
    tiempo_intervalo = 111111 #Como si me dijeran que de 0 a 111111 segundos entra la masa x (lo que entrgara el return)
    if t < tiempo_intervalo:
        return M_IN 
    else:
        return 0.0 #en caso de que no hay entrada debe retornar 0 

def flujo_salida_func(t): #aca hay que hacer lo mismo que arriba
    return M_OUT #por ahora lo dejo constante




#Ahora para correr la simulacion

#Tiempo que quiero que este simulando 
T_SIMULACION = 1111  #(s)
#Tiempo para la simulacion como rango 
t_span = [0, T_SIMULACION]

#Para guardar los resultados y despues graficar 
t_eval = np.linspace(t_span[0], t_span[1], 500)


#Todos los parametros juntos para pasarlos al solver
args_modelo = (V, AREA_S, U_TW, T_AMB, H_IN,M_IN, M_OUT, FLUIDO)


#Ahora llamamos al solver de SciPy
solucion = solve_ivp(modelo_tanque, t_span,y0,                    
    method='BDF',           #Lo recomendaban para sistemas stiff (rigidos)
    t_eval=t_eval,args=args_modelo)



#Para obtener los resultados y graficar 
if solucion.success:
    t = solucion.t
    T = solucion.y[0]  
    P_pascales = solucion.y[1] 
    P_bar = P_pascales / 1e5   #para pasar de Pa a bar para graficar
    

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    #Temperatura vs Tiempo
    ax1.plot(t, T, label='Temperatura (K)', color='red')
    ax1.set_ylabel('Temperatura (K)')
    ax1.legend()
    ax1.grid(True)
    ax1.set_title(f'Simulación Tanque de Hidrogeno')
    
    #Presión vs Tiempo
    ax2.plot(t, P_bar, label='Presión (bar)', color='blue')
    ax2.set_xlabel('Tiempo (s)')
    ax2.set_ylabel('Presión (bar)')
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    plt.show()

else:
    print("La simulación falló")
    print(solucion.message)