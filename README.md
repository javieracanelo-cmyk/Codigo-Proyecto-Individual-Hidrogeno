En el repositorio se entrega un modelo termodinámico programado en python para simular el acoplado de un tanque de hidrógeno
comprimido y un intercambiador de calor de placas corrugadas 1D a lo largo de la longitud y en función del tiempo con configuración contracorriente. 

El objetivo era estudiar el comportamiento transiente del sistema de almacenamiento y recuperación de calor para una planta de rSOC. Se buscaba lograr una simulación detallada y modificable para evaluar la implementación de un intercambiador de calor en reemplazo de un calefactor eléctrico y con uso de una corriente de deshecho de aire previamente calentada. 
El principio de operación se dividirá en dos partes: 

##Tanque de H₂

El tanque se modela como un sistema transiente de volumen rígido, perfectamente mezclado, que opera solamente con salida de flujo constante. La dinámica se describe mediante balances diferenciales de masa y energía directamente de la formulación propuesta por Reznicek y Braun (2018). 
Se utiliza un base molar, y se resuelve obteniendo las propiedades del fluido con la libreria CoolProp y simulando la presión con la ecuación de estado cúbica Redlich-Kwong (RK). 
Ecuación de Redlich-Kwong: 

$$P = \frac{R T}{v-b} - \frac{a}{\sqrt{T} v(v+b)}$$

donde, 

$$
v = \frac{V}{n}
$$

y

$$
a = 0.42748\ \frac{R^{2}\ T_c^{2.5}}{P_c}, \qquad 
b = 0.08664\\frac{R\ T_c}{P_c}.
$$

El balance de masa se define como: 

$$\frac{dn}{dt}=-\dot{n}_{out}$$

El balance de energía se define como: 

$$\frac{dU}{dt}=\dot{H}_{in}-\dot{H}_{out}-\dot{Q}_{loss}$$

donde, 

$$\dot{Q}_{loss}=U_{tw}\,A_s\,(T-T_{amb})$$

con: 

$$U_{tw} = 5 \frac{W}{m^2K}$$

$$T_{amb} = 25 °C$$ 

Para el Área superficial del tanque se considera en área total de un cilindro: 

$$2\pi r L + 2 \pi r^2$$ 

Considerando la relación para el largo a partir del diametro: 

$$L = 3D$$

Obtenida en Reznicek y Braun (2018). Para el cálculo de la masa incial se ocupó la fórmula: 

$$I_{\mathrm{tot}}=i_{\mathrm{fuelcell}}A_{\mathrm{cell}}N_{\mathrm{cell}}N_{\mathrm{stacks}}N_{\mathrm{arrays}}N_{\mathrm{mod}}$$

$$n_0=\frac{\lvert v_{\mathrm{H}_2}\rvert I_{\mathrm{tot}}\Delta t}{n_S S F U_F}$$

$$m_0=n_0 M_{\mathrm{H}_2}$$

##Intercambiador de Calor 
Para el caso del intercambiador de calor, se utiliza un intercambiador de tipo de placas corrugadas por la necesidad de un aparato compacto para la planta de rSOC. 
El intercambiador se modela como un sistema 1D en el eje longitudinal, en N segmentos, para tener consistencia física. En cada segmento se resuelven balances difreneciales de energía para el fluido frío (hidrógeno), fluido claiente (aire) y la pared metálica. 

Los balances por segmento

$$\frac{dT_c}{dt} = \frac{\dot{m}_c c_{p,c}(T_c^{\text{in}} - T_c^{\text{out}}) + A_{\text{transfer}} h_{\mathrm{conv},c}(T_w - T_c)}{m_c c_{p,c}}$$

$$\frac{dT_h}{dt} = \frac{\dot{m}_h c_{p,h}(T_h^{\text{in}} - T_h^{\text{out}}) - A_{\text{transfer}} h_{\mathrm{conv},h}(T_h - T_w)}{m_h c_{p,h}}$$

$$\frac{dT_w}{dt} = \frac{A_{\text{transfer}} \left[ h_{\mathrm{conv},h}(T_h - T_w) - h_{\mathrm{conv},c}(T_w - T_c) \right]}{C_w}$$

Se define el esquema como contracorriente, para tener una mayor eficiencia en la transferencia de calor entre los fluidos. De esa forma la entrada del hidrógeno se define en x = 0 y la entrada del aire en x = L (recorriendo ambos fluidos los segmentos en sentido inverso). 
Para calcular el área de transferencia se establecen las dimensiones físicas, los materiales y la configuración de canales. 

Propiedades Físicas del Material

| Parámetro | Símbolo | Valor | Unidad | Descripción | 
| :--- | :--- | :--- | :--- | :--- |
| Calor específico del metal | $c_{p,w}$ | $500$ | $\frac{J}{kg K}$ | Propiedad del acero inoxidable de las placas. |
| Densidad del metal | $\rho_{w}$ | $8000$ | $\frac{kg}{m^3}$ | Densidad del acero. |
| Espesor de las placas | $e$ | $0.0008$ | m | Espesor de cada placa. |

Configuración Geométrica y Canales

| Parámetro | Símbolo | Valor | Unidad | Descripción |
| :--- | :--- | :--- | :--- | :--- |
| Número total de placas | $N_{placas}$ | $180$ | - | Número total de placas ( se ocupa el máximo). |
| Longitud de las placas | $L$ | $0.381$ | m | Longitud del canal de flujo. |
| Ancho de las placas | $W$ | $0.175$ | m | Ancho efectivo del canal. |
| Distancia entre chevrones | $p$ | $0.01$ | m | Distancia entre corrugaciones. |
| Ángulo chevron | $\beta$ | $10.041^{\circ}$ | rad | Ángulo de las corrugaciones. |
| Espacio entre placas | $b$ | $0.0014352$ | m | Espacio libre para el fluido ($0.0022352 - e$). |

Los parámetros se obtuvieron a partir de la selección de geometría aplicada anteriormente en intercambiadores de este tipo  (Al-Dahhan et al., 2024), consistente con las condiciones de la simulación presente.

El número total de canales es $N_{canales} = N_{placas} - 1$, el cual se divide en $N_{frio}$ (canales de H₂) y $N_{caliente}$ (canales de Aire) de forma alternada.

Cálculo de Parámetros Hidráulicos y de Transferencia

1. Factor de Ampliación por Corrugación ($\phi$)

El factor de corrugación ($\phi$) representa la relación entre el área superficial real y el área superficial plana, debido a la geometría de las placas corrugadas. Se calcula a partir de la razón de corrugación $\gamma$ ($\gamma = 2b/p$) mediante una correlación empírica:

$$\phi = \frac{1}{6} \left( 1 + \sqrt{1 + \gamma^2 \left(\frac{\pi}{2\cos\beta}\right)^2} + 4 \sqrt{1 + \gamma^2 \left(\frac{\pi}{2\sqrt{2}\cos\beta}\right)^2} \right)$$

2. Diámetro Hidráulico Corregido ($\mathbf{D_h}$)

El diámetro hidráulico corregido considera el efecto de la corrugación en el flujo:

$$D_h = \frac{2b}{\phi}$$

3. Áreas de Flujo y Volumen

El área transversal total ($A_{cs}$) para cada fluido y el volumen total ($V$) del fluido dentro del ICC son:

$$A_{cs, c} = (W \cdot b) \cdot N_{frio}$$

$$V_c = A_{cs, c} \cdot L$$

4. Área Real de Transferencia ($\mathbf{A_{real}}$)

El área total disponible para la transferencia de calor entre fluidos, ajustada por el factor $\phi$ y excluyendo las dos placas de los extremos:

$$A_{\text{real}} = (N_{placas} - 2) \cdot W \cdot L \cdot \phi$$

5. Capacidad Térmica de la Pared ($\mathbf{C_w}$)

La capacidad térmica total de la masa de metal del intercambiador se calcula como:

$$C_w = m_{tot,w} \cdot c_{p,w} \quad \text{donde} \quad m_{tot,w} = \rho_w \cdot e \cdot W \cdot L \cdot N_{placas} \cdot \phi$$


Presión de Operación 

Para el cálculo de la presión se trabajó a presión constante considerando una presión promedio entre la presión incial y la presión final, donde se consideró una caida de presión de 500 Pa a partir de lo dicho por Zang et al. (2023), en el análisis de una planta Power-to-X flexible basada en celdas de óxido sólido reversibles. 


Cálculo de Coeficientes de Transferencia de Calor por Convección ($\mathbf{h_{conv}}$)

Para simular la transferencia de calor entre el fluido y la pared, es necesario calcular el coeficiente de transferencia de calor por convección ($h_{conv}$) en cada canal. Este cálculo depende del **régimen de flujo** (Laminar o Turbulento), definido por el número de Reynolds.

1. Propiedades Físicas y Adimensionales (H₂ y Aire)

Ambas funciones obtienen las propiedades del fluido (densidad $\rho$, viscosidad $\mu$, calor específico $c_p$ y conductividad térmica $k$) en las condiciones de presión ($P$) y temperatura ($T$) actuales, utilizando la librería CoolProp.

Con estas propiedades, se calculan los números adimensionales clave:

* **Velocidad Promedio ($\mathbf{u}$):** $u = \dot{m} / (\rho \cdot A_{cs})$
* **Número de Reynolds ($\mathbf{Re}$):** Define el régimen de flujo.
    $$\mathrm{Re} = \frac{\rho \cdot u \cdot D_h}{\mu}$$
* **Número de Prandtl ($\mathbf{Pr}$):** Relaciona la difusión de momento con la difusión térmica.
    $$\mathrm{Pr} = \frac{c_p \cdot \mu}{k}$$

2. Correlaciones del Número de Nusselt ($\mathbf{Nu}$)

El coeficiente $h_{conv}$ se obtiene a través del **Número de Nusselt ($\mathrm{Nu}$)**, para el cual se aplican diferentes correlaciones según el fluido:

A. Régimen Laminar ($\mathbf{Re < 2300}$)

Para el flujo laminar, se utiliza una correlación de Nusselt adecuada para canales de intercambiadores de placas que tienen en cuenta el desarrollo de la longitud ($L_{adim}$) anteriormente mencionada por Naphon et al. (2018):

$$\mathrm{Nu} = 3.63 + 0.086 \cdot \frac{L_{\text{adim}}^{-1.33}}{1 + 0.1 \cdot \mathrm{Pr} \cdot (\frac{D_h \cdot \mathrm{Re}}{L})^{0.83}}$$

Donde $L_{adim}$ es la longitud adimensional: $L_{\text{adim}} = \frac{L/D_h}{\mathrm{Re} \cdot \mathrm{Pr}}$



B. Régimen Turbulento ($\mathbf{Re \geq 2300}$)

Para el flujo turbulento, se emplea la correlación de Gnielinski, donde el factor de fricción ($f$) es proporcionado por la ecuación de Konakov (Konakov, 1946):

$$f = 0.25 \cdot (1.82 \cdot \ln(\mathrm{Re}) - 1.5)^{-2}$$

El Número de Nusselt se calcula como:

$$\mathrm{Nu} = \frac{(f/2) \cdot (\mathrm{Re} - 1000) \cdot \mathrm{Pr}}{1 + 12.7 \sqrt{f/2} \cdot (\mathrm{Pr}^{2/3} - 1)}$$

3. Cálculo Final Coeficiente de Transferencia de Calor ($\mathbf{h_{conv}}$)

Finalmente, el coeficiente de transferencia de calor por convección se calcula a partir del Número de Nusselt y la conductividad térmica del fluido:

$$h_{\mathrm{conv}} = \frac{\mathrm{Nu} \cdot k}{D_h}$$



El código obtenido entrega una gráfica completa de la variación de temperatura, presión y masa dentro del estanque con respecto al tiempo. Junto a una imagen completa de la variación de temperatura para cada fluido a lo largo del intercambiador. 

Finalmente, como se dice en las noticias actuales Chile se posiciona como un líder mundial en la producción de Hidrógeno Verde, sobre todo por su potencial de energías renovables. La optimización de la eficiencia en las plantas de rSOC, hace más atractiva la inversión económica en estas tecnologías.
La versatilidad de estas plantas, permite su uso tanto en la electrolisis, en el almacenamiento de energías renovables (frente a la alta intermitencia de la energía solar y eólica) y como suministro de energía eléctrica al hacer uso del hidrógeno en celda de combustible. 
Su factibilidad depende del avance de la industria del Hidrógeno verde en el país.




Referencias

Al-Dahhan, Z. S. H., Al-Shaghdali, H. M., Al-Jubouri, M. H., Al-Shukri, M. M., A. M. H. Al-Dulaimi, N. A. H. Al-Anssari, & Al-Saadi, B. A. M. (2024). Thermal and hydrodynamic optimization of plate heat exchanger using multi objective genetic algorithm. International Communications in Heat and Mass Transfer, 150, 107227.

Konakov, P. K. (1946). Eine neue Formel für den Reibungskoeffizienten glatter Rohre (Una nueva fórmula para el coeficiente de fricción de tubos lisos). Berichte Akad, Wiss. UdSSR., 51, 503–506. (Obra original publicada en ruso).

Naphon, M. R., Wongwises, S., da Silva, A. L. N., & S. Z. N. O. (2018). Thermal and hydrodynamic analysis of a cross-flow compact heat exchanger. International Journal of Heat and Mass Transfer, 117, 300–311.

Reznicek, E., & Braun, R. J. (2018). Techno-economic and off-design analysis of stand-alone, distributed-scale reversible solid oxide cell energy storage systems. Energy Conversion and Management, 175, 263–277. Elsevier.

Zang, L., Zhou, G., Luo, X., Wang, S., Lee, Y. J., da Silva, C. A. G. O., de Oliveira, R. A. F., P. A. N. S. Filho, D. B. (2023). Energy and environmental analysis of a flexible Power-to-X plant based on Reversible Solid Oxide Cells (rSOCs) for an urban district. Energy, 278, 127929.


