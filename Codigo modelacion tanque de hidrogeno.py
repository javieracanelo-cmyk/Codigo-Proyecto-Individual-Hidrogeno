#Codigo Tanque 
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

#El modelo del tanque se define a nivel molar y con respecto a T, ocupando coolprop 

#Para implementar RK para P 
# Parámetros Redlich-Kwong para H2
R = 8.314  # J/(mol·K)  

Tc_H2 = 33.19      # K  (temperatura crítica del H2, aprox)
Pc_H2 = 1.296e6    # Pa (presión crítica del H2, aprox)

# Parámetros RK
a_RK = 0.42748 * R**2 * Tc_H2**2.5 / Pc_H2
b_RK = 0.08664 * R * Tc_H2 / Pc_H2

def P_RK(T, n, V):
    """
    Presión Redlich-Kwong:
    P = R T / (v - b) - a / (sqrt(T) v (v + b)),  v = V/n
    """
    v = V / n  # volumen molar [m³/mol]
    return R*T/(v - b_RK) - a_RK/(np.sqrt(T)*v*(v + b_RK))
#Función de modelación: 
def modelo_tanque_py(t, y, V, Utw, T_amb, M_H2, R, T_in, P_in, m_dot_in, m_dot_out, T0):
    
    #Funcion del modelo del tanque ocupando: 
    #y(1) = n
    #y(2) = n*h
    #h=y(2)/y(1)
    
    #Parametros del tanque:
    d= (4 * V / (3 * np.pi))**(1/3) 
    r = d / 2
    L = 3 * d #A partir del articulo 
    As = 2*np.pi*r**2 + 2*np.pi*r*L #Área superficial

    #Definir el vector 
    n = y[0]     #moles
    nh = y[1]   #Entalpía (J)
    
    if n < 1e-6: #como ocupamos una división por n, esto evita que se de una división por 0 
        n = 1e-6 
        nh = 0

    #Definicion del sistema como queremos (en molar)
    h = nh / n        #Entalpía molar (J/mol) -> h = y(2)/y(1)
    h_mass = h / M_H2  #Entalpía másica (J/kg) (así se ocupa CoolProp)
    v_m = V / n        #Volumen molar (m³/mol)

    #Entalpías calculadas con Coolprop (esto es para cuando agreguemos una posible entrada y para definir la salida)
    h_in_mass = PropsSI('Hmass', 'T', T_in, 'P', P_in, 'Hydrogen') #J/kg
    h_in = h_in_mass * M_H2 #J/mol (conversión a molar)
    h_out = h #J/mol (Se define que la corriente de salida tiene las mismas propiedades que el tanque)

    #Iteración para los calculos de T y P -> Explicación del código: Yo le doy T0 y P0 con eso calcula n0 y h0
    #El solver ira resolviendo dn/dt y dn*h/dt con respecto a los balances y con n y h en cada momento calcula T y P
    Titer = T0  #supuesto inicial para comenzar la iteración 
    tol = 1e-3
    P_tanque = 0 #para crear la variable de la presión 
    
    for _ in range(20):
        
        #La presión se calcula con Redlich-Kwong
        P_tanque = P_RK(Titer, n, V)
        
        
        try:
            #Aca ocupo Coolprop para obtener la temperatura a partir de h y P 
            Tnew = PropsSI('T', 'Hmass', h_mass, 'P', P_tanque, 'Hydrogen')
        except ValueError: #En caso de que Coolprop falle, se mantiene la temperatura anterior (Acá debería ir lo del NIST pero no logro aplicarlo)
            Tnew = Titer
            break

        if abs(Tnew - Titer) < tol: #Si las dos T son iguales es correcto y entra al siguiente paso
            break
            
        Titer = Tnew
        
    T_tanque = Titer 
    P_tanque = R * T_tanque / v_m

    #Aqui se definen los balances de masa y energía
    dn_dt = (m_dot_in / M_H2) - (m_dot_out / M_H2)
    Q_loss = Utw * As * (T_tanque - T_amb)
    
    dn_h_dt = (m_dot_in / M_H2 * h_in) - \
              (m_dot_out / M_H2 * h_out) - \
              Q_loss

    #Para que me devuelva las derivadas 
    return [dn_dt, dn_h_dt]

#Ahora para simular el tanque en la planta rSOC 

#Calculo de la masa (para que no me quede sin hidrógeno en el tanque en el tiempo de funcionamiento y según el requerimiento de la celda)
M_H2 = 2.016e-3   #kg/mol (Masa molar hidrógeno)
Deltat = 4*3600 #s (Tiempo en el que quiero que esté operando como celda de combustible)
v_h2 = -1 #coeficiente estequeométrico (negativo porque es consumo de H2)
i_fuelcell = 3500 #A/m² Corriente de la celda de combustible
A_cell = 0.0088208 #m² Área activa de la celda de combustible (donde se aplica la corriente)
n = 2 #n° de electrones transferidos por mol de H2
UF = 0.7 #Factor de utilización del hidrógeno
Ncell = 68 #Número de celdas por stack
Nstacks = 10 #Número de stacks por arreglo 
Narrays = 10 #Número de arreglos por módulo
Nmod = 1 #Número de módulos

m0 = Deltat*((abs(v_h2)*i_fuelcell*A_cell*M_H2)/(n*96485*UF))*Ncell*Nstacks*Narrays*Nmod
n0 = m0 / M_H2 #kg de hidrógeno inicial en el tanque

#Condiciones iniciales 
T0 = 273.15 + 40      #K
P0 = 70e5         #Pa

#Constantes del tanque 
rho0 = PropsSI('D', 'T', T0, 'P', P0, 'Hydrogen') #kg/m³ (Densidad a esa T0 Y P0)
V = m0/rho0          #m³ (Calculado con la ecuación)
Utw = 5.0           #W/(m²·K) (A partir del articulo)
T_amb = 25 + 273.15 #K (Condiciones estándar)



#En caso de que haya una entrada acá se debería agregar (para el futuro)
T_in = 293.15     #K 
P_in = 1e5        #Pa 
m_dot_in = 0.0    #kg/s  (Así se define sin entrada)
m_dot_out = 0.005  #kg/s  (Hay salida pequeña y constante)

#Calculo de entalpía inicial
h0_mass = PropsSI('Hmass', 'T', T0, 'P', P0, 'Hydrogen') #J/kg (entalpía inicial con Coolprop)
h0_mol = h0_mass * M_H2 #J/mol (Conversión a molar)
n_h0 = n0 * h0_mol      #J

#Juntamos los moles y la entalpía inicial en el vector y0
y0 = [n0, n_h0] 

#Inicio de la modelación
t_span = [0, 7200] #Para definir por cuanto tiempo se simule
t_eval = np.linspace(t_span[0], t_span[1], 300) #Para ir viendo T y P cada segundo

args_sim = (V, Utw, T_amb, M_H2, R, T_in, P_in, m_dot_in, m_dot_out,T0)

sol = solve_ivp(
    modelo_tanque_py,
    t_span,
    y0,
    method='BDF',
    t_eval=t_eval,
    args=args_sim
)


#Para graficar T y P con respecto al tiempo
#El solver retorna las derivadas del balance de masa, por lo que necesitamos calcular a partir de eso P y T
t = sol.t #tiempos
n_t = sol.y[0] #extrae los moles en cada punto
n_h_t = sol.y[1] #extrae las entalpías en cada punto

T_t = np.zeros_like(n_t) #Arreglo para guardar las temperaturas
P_t = np.zeros_like(n_t) #Arreglo para guardar las presiones
T_resp = T0 #Temperatura en t = 0

for i in range(len(t)): #Para sacar T y P en cada tiempo
    n_i = n_t[i]
    n_h_i = n_h_t[i]
    
    #Evitar división por cero si el tanque está vacío
    if n_i < 1e-6:
        T_t[i] = T_amb
        P_t[i] = 0.0
        continue
        
    h_mol_i = n_h_i / n_i 
    v_m_i = V / n_i
    h_mass_i = h_mol_i / M_H2
    
    #Para encontrar T y P a partir de n y h, se itera igual que en el modelo
    Titer = T_resp
    for _ in range(20):
         #Presión con Redlich-Kwong
        Piter = P_RK(Titer, n_i, V)
        try:
            Tnew = PropsSI('T', 'Hmass', h_mass_i, 'P', Piter, 'Hydrogen')
        except ValueError:
            Tnew = Titer
            break
        
        if abs(Tnew - Titer) < 1e-3:
            break
        Titer = Tnew
        
    T_t[i] = Titer
    P_t[i] = (n_i * R * Titer) / V
    T_resp = Titer #Actualiza la T

m_t = n_t * M_H2  #Masa en kg

#Gráficos: T(t), P(t), m(t)

# 1) Temperatura 
plt.figure(figsize=(8,4))
plt.plot(t, T_t, color='red', linewidth=2, label='T (K)')
plt.xlabel('Tiempo (s)'); plt.ylabel('Temperatura (K)')
plt.title('Temperatura del tanque'); plt.grid(True); plt.legend()
plt.tight_layout(); plt.show()

# 2) Presión 
plt.figure(figsize=(8,4))
plt.plot(t, P_t/1e5, color='#1f3b77', linestyle='--', linewidth=2, label='P (bar)')
plt.xlabel('Tiempo (s)'); plt.ylabel('Presión (bar)')
plt.title('Presión del tanque'); plt.grid(True); plt.legend()
plt.tight_layout(); plt.show()

# 3) Masa 
plt.figure(figsize=(8,4))
plt.plot(t, m_t, color='green', marker='o', markevery=20, linewidth=1.8, label='m (kg)')
plt.xlabel('Tiempo (s)'); plt.ylabel('Masa H$_2$ (kg)')
plt.title('Masa de H$_2$ en el tanque'); plt.grid(True); plt.legend()
plt.tight_layout(); plt.show()
