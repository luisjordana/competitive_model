# Competitive Model

In this project, modeling outcomes following Kelso et al. () will be presented to exemplify how antagonistic coupling affects rhythmical coupled behavior. This specific way of modeling follows the formula presented below:

$$x ̈+(αx^2+βx ̇^2-γ) x ̇+ω^2 x=(A+B(x+μ_1 y)^2)(x ̇-μ_1 y ̇)$$

$$y ̈+(αy^2+βy ̇^2-γ) y ̇+ω^2 y=(A+B(y+μ_2 x)^2)(y ̇-μ_2 x ̇)$$

Such that x, ̇x and ̈x represent the position, velocity, and acceleration, respectively, of oscillator one, while the y, ̇y and ̈y represent the position and velocity of oscillator two. On the other hand, α, β, γ, ω^2, A and B represent constants. Constants α, and β, control the nonlinear damping of the individual oscillators, while γ controls the linear damping of such oscillators. Furthermore, ω^2 controls the period of oscillation. Finally, the values of A and B, as well as the relation between them, control the stability of the system. In symmetrical attractive settings (when the behavior of the two oscillators attract each other) this means that changing the relation between A and B controls if the system is in a bi-stable (0 º and 180 º) or a mono-stable (0 º) regime. Following Kelso et al.13, the novelty of this modeling comes from parameter μ. Such a parameter, when set to a negative value in one oscillator, allows generating a repulsive force, which makes for a repulsive-attractive coupling, which can be used to model antagonistic coupling.

The project was done with Matlab and includes a Simulink model (image of the model can be seen in Simulink_Model.png), images of different coupling outcomes as well as some necessary functions, and  run_hkb_coupled_osc_fuchs that was used to run the Simulink model with different values of each constant leading to different coupling results.

