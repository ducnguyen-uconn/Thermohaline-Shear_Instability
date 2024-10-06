Governing equations ([Radko, 2016](https://doi.org/10.1017/jfm.2016.547)) for perturbation of velocity $\boldsymbol{u}$, temperature $T$, and salinity $S$:

$$\nabla\cdot \textbf{u} = 0$$

$$\partial_t \boldsymbol{u} + U_{bg}\partial_x\boldsymbol{u} + w\partial_z U_{bg}\boldsymbol{e}_x + \boldsymbol{u}\cdot\nabla\boldsymbol{u} = -\nabla p + \frac{Pr}{Pe}\nabla^2\boldsymbol{u} + \frac{4\pi^2 Ri}{R_\rho-1}(T-S)\boldsymbol{e}_z,$$

$$\partial_t T + U_{bg}\partial_x T + \boldsymbol{u}\cdot\nabla T - w = \frac{1}{Pe}\nabla^2 T,$$

$$\partial_t S + U_{bg}\partial_x S + \boldsymbol{u}\cdot\nabla S - R_\rho w = \frac{\tau}{Pe}\nabla^2 S$$

where $U_{bg}=sin(2\pi z)$.

For example, with 3D configuration, these governing equations can be expressed as follows

$$\nabla\cdot \boldsymbol{u} = 0,$$

$$\partial_t u + sin(2\pi z)\partial_x u +2\pi cos(2\pi z)w+ \boldsymbol{u}\cdot\nabla u = -\partial_x p + \frac{Pr}{Pe}\nabla^2u ,$$

$$\partial_t v + sin(2\pi z)\partial_x v + \boldsymbol{u}\cdot\nabla v = -\partial_y p + \frac{Pr}{Pe}\nabla^2v ,$$

$$\partial_t w + sin(2\pi z)\partial_x w + \boldsymbol{u}\cdot\nabla w = -\partial_z p + \frac{Pr}{Pe}\nabla^2w + \frac{4\pi^2 Ri}{R_\rho-1}(T-S),$$

$$\partial_t T + sin(2\pi z)\partial_x T + \boldsymbol{u}\cdot\nabla T - w = \frac{1}{Pe}\nabla^2 T,$$

$$\partial_t S + sin(2\pi z)\partial_x S + \boldsymbol{u}\cdot\nabla S - R_\rho w = \frac{\tau}{Pe}\nabla^2 S.$$

Now, we decompose perturbations of velocity, pressure, temperature, and salinity fields into two components: mean and fluctuation quantities.

$$u = \bar{u}+u'$$

$$w = \bar{w}+w' = w', (\bar{w}=0)$$

$$p = \bar{p}+p',$$

$T = \bar{T}+T',$

$S = \bar{S}+S'.$

Then, we have

$$
\boldsymbol{\bar{u}}\cdot\nabla\boldsymbol{\bar{u}} = 0,\\
<\boldsymbol{u}'\cdot\nabla\boldsymbol{\bar{u}}>_h = <\boldsymbol{u}'>_h\cdot\nabla\boldsymbol{\bar{u}} = 0,\\
<\boldsymbol{\bar{u}}\cdot\nabla\boldsymbol{u}'>_h = \boldsymbol{\bar{u}}\cdot\nabla<\boldsymbol{u}'>_h = 0,\\
<\boldsymbol{u}'\cdot\nabla\boldsymbol{u}'>_h  \neq 0 ,
$$

and 
$$
\boldsymbol{\bar{u}}\cdot\nabla\boldsymbol{\bar{T}} = 0,\\
<\boldsymbol{u}'\cdot\nabla\boldsymbol{\bar{T}}>_h = <\boldsymbol{u}'>_h\cdot\nabla\boldsymbol{\bar{T}} = 0,\\
<\boldsymbol{\bar{u}}\cdot\nabla\boldsymbol{T}'>_h = \boldsymbol{\bar{u}}\cdot\nabla<\boldsymbol{T}'>_h = 0,\\
<\boldsymbol{u}'\cdot\nabla\boldsymbol{T}'>_h  \neq 0 ,
$$

and 
$$
\boldsymbol{\bar{u}}\cdot\nabla\boldsymbol{\bar{S}} = 0,\\
<\boldsymbol{u}'\cdot\nabla\boldsymbol{\bar{S}}>_h = <\boldsymbol{u}'>_h\cdot\nabla\boldsymbol{\bar{S}} = 0,\\
<\boldsymbol{\bar{u}}\cdot\nabla\boldsymbol{S}'>_h = \boldsymbol{\bar{u}}\cdot\nabla<\boldsymbol{S}'>_h = 0,\\
<\boldsymbol{u}'\cdot\nabla\boldsymbol{S}'>_h  \neq 0 ,
$$