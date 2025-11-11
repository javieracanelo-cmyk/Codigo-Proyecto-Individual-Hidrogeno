#Codigo Tanque 
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

#El modelo del tanque se define a nivel molar y con respecto a T, ocupando coolprop y asumiendo gas ideal para la presión
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
    As = 2*np.pi*r**2 + 2*np.pi*r*L # Área transversal 

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
        
        #La presión se calcula como gas ideal
        P_tanque = R * Titer / v_m
        
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

#Constantes del tanque 
V = 14.3          #m³ (Calculado con la ecuación)
Utw = 5.0           #W/(m²·K) (A partir del articulo)
T_amb = 25 + 273.15 #K (Condiciones estándar)
M_H2 = 2.016e-3   #kg/mol (Masa molar hidrógeno)
R = 8.314       #J/(mol·K)

#En caso de que haya una entrada acá se debería agregar (para el futuro)
T_in = 293.15     #K 
P_in = 1e5        #Pa 
m_dot_in = 0.0    #kg/s  (Así se define sin entrada)
m_dot_out = 0.01  #kg/s  (Hay salida pequeña y constante)

#Condiciones iniciales 
T0 = 273.15 + 40      #K
P0 = 70e5         #Pa

#Calculo de moles y entalpía inicial
rho0 = PropsSI('D', 'T', T0, 'P', P0, 'Hydrogen') #kg/m³ (Densidad a esa T0 Y P0)
n0 = (rho0 * V) / M_H2  #moles a partir de la densidad y el volumen 
h0_mass = PropsSI('Hmass', 'T', T0, 'P', P0, 'Hydrogen') #J/kg (entalpía inicial con Coolprop)
h0_mol = h0_mass * M_H2 #J/mol (Conversión a molar)
n_h0 = n0 * h0_mol      #J

#Juntamos los moles y la entalpía inicial en el vector y0
y0 = [n0, n_h0] 

#Inicio de la modelación
t_span = [0, 300] #Para definir por cuanto tiempo se simule
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
        Piter = (n_i * R * Titer) / V # Gas Ideal
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

#Gráfico
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax1.plot(t, T_t, label='Temperatura (K)', color='red')
ax1.set_ylabel('Temperatura (K)')
ax1.legend()
ax1.grid(True)
ax1.set_title(f'Simulación Vaciado de Tanque (Gas Ideal + CoolProp H)')

ax2.plot(t, P_t / 1e5, label='Presión (bar)', color='blue') # Convertir Pa a bar
ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel('Presión (bar)')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()