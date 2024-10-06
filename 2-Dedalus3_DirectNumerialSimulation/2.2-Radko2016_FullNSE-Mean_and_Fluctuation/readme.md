Governing equations ([Radko, 2016](https://doi.org/10.1017/jfm.2016.547)) for perturbation of velocity $\vec{u}$, temperature $T$, and salinity $S$:

$$\nabla\cdot \vec{u} = 0$$

$$\partial_t \vec{u} + U_{bg}\partial_x\vec{u} + w\partial_z U_{bg}\vec{e}_x + \vec{u}\cdot\nabla\vec{u} = -\nabla p + \frac{Pr}{Pe}\nabla^2\vec{u} + \frac{4\pi^2 Ri}{R_\rho-1}(T-S)\vec{e}_z,$$

$$\partial_t T + U_{bg}\partial_x T + \vec{u}\cdot\nabla T - w = \frac{1}{Pe}\nabla^2 T,$$

$$\partial_t S + U_{bg}\partial_x S + \vec{u}\cdot\nabla S - R_\rho w = \frac{\tau}{Pe}\nabla^2 S$$

where $U_{bg}=sin(2\pi z)$.

For example, with 3D configuration, these governing equations can be expressed as follows

$$\nabla\cdot \vec{u} = 0,$$

$$\partial_t u + sin(2\pi z)\partial_x u +2\pi cos(2\pi z)w+ \vec{u}\cdot\nabla u = -\partial_x p + \frac{Pr}{Pe}\nabla^2u ,$$

$$\partial_t v + sin(2\pi z)\partial_x v + \vec{u}\cdot\nabla v = -\partial_y p + \frac{Pr}{Pe}\nabla^2v ,$$

$$\partial_t w + sin(2\pi z)\partial_x w + \vec{u}\cdot\nabla w = -\partial_z p + \frac{Pr}{Pe}\nabla^2w + \frac{4\pi^2 Ri}{R_\rho-1}(T-S),$$

$$\partial_t T + sin(2\pi z)\partial_x T + \vec{u}\cdot\nabla T - w = \frac{1}{Pe}\nabla^2 T,$$

$$\partial_t S + sin(2\pi z)\partial_x S + \vec{u}\cdot\nabla S - R_\rho w = \frac{\tau}{Pe}\nabla^2 S.$$

Now, we decompose perturbations of velocity, pressure, temperature, and salinity fields into two components: mean and fluctuation quantities.

$$u = \bar{u}+u'$$

$$w = \bar{w}+w' = w', (\bar{w}=0)$$

$$p = \bar{p}+p',$$

$T = \bar{T}+T',$

$S = \bar{S}+S'.$

Then, we have

$$
\vec{\bar{u}}\cdot\nabla\vec{\bar{u}} = 0,\\
<\vec{u}'\cdot\nabla\vec{\bar{u}}>_h = <\vec{u}'>_h\cdot\nabla\vec{\bar{u}} = 0,\\
<\vec{\bar{u}}\cdot\nabla\vec{u}'>_h = \vec{\bar{u}}\cdot\nabla<\vec{u}'>_h = 0,\\
<\vec{u}'\cdot\nabla\vec{u}'>_h  \neq 0 ,
$$

and 
$$
\vec{\bar{u}}\cdot\nabla\vec{\bar{T}} = 0,\\
<\vec{u}'\cdot\nabla\vec{\bar{T}}>_h = <\vec{u}'>_h\cdot\nabla\vec{\bar{T}} = 0,\\
<\vec{\bar{u}}\cdot\nabla\vec{T}'>_h = \vec{\bar{u}}\cdot\nabla<\vec{T}'>_h = 0,\\
<\vec{u}'\cdot\nabla\vec{T}'>_h  \neq 0 ,
$$

and 
$$
\vec{\bar{u}}\cdot\nabla\vec{\bar{S}} = 0,\\
<\vec{u}'\cdot\nabla\vec{\bar{S}}>_h = <\vec{u}'>_h\cdot\nabla\vec{\bar{S}} = 0,\\
<\vec{\bar{u}}\cdot\nabla\vec{S}'>_h = \vec{\bar{u}}\cdot\nabla<\vec{S}'>_h = 0,\\
<\vec{u}'\cdot\nabla\vec{S}'>_h  \neq 0 ,
$$