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

$$v = \bar{v}+v'$$

$$w = \bar{w}+w' = w', (\bar{w}=0)$$

$$p = \bar{p}+p',$$

$$T = \bar{T}+T',$$

$$S = \bar{S}+S'.$$

Based on this, we have horizontally averaged quantities

$$
\langle\vec{u}\rangle_h=\bar{\vec{u}}(z), \quad \langle\bar{T}\rangle_h=\bar{T}(z), \quad \langle\bar{S}\rangle_h=\bar{S}(z), \\
\langle\vec{u}'\rangle_h=0, \quad \langle T'\rangle_h=0, \quad \langle S'\rangle_h=0,
$$

where $\langle.\rangle_h = \frac{1}{L_x L_y}\int\langle.\rangle dxdy$

Then, we have

$$
\vec{\bar{u}}\cdot\nabla\vec{\bar{u}} = 0,\\
\Rightarrow\langle\vec{u}'\cdot\nabla\vec{\bar{u}}\rangle_h = \langle\vec{u}'\rangle_h\cdot\nabla\vec{\bar{u}} = 0,\\
\Rightarrow\langle\vec{\bar{u}}\cdot\nabla\vec{u}'\rangle_h = \vec{\bar{u}}\cdot\nabla\langle\vec{u}'\rangle_h = 0,\\
\Rightarrow\langle\vec{u}'\cdot\nabla\vec{u}'\rangle_h  \neq 0 ,
$$

and 

$$
\vec{\bar{u}}\cdot\nabla\vec{\bar{T}} = 0,\\
\Rightarrow\langle\vec{u}'\cdot\nabla\vec{\bar{T}}\rangle_h = \langle\vec{u}'\rangle_h\cdot\nabla\vec{\bar{T}} = 0,\\
\Rightarrow\langle\vec{\bar{u}}\cdot\nabla\vec{T}'\rangle_h = \vec{\bar{u}}\cdot\nabla\langle\vec{T}'\rangle_h = 0,\\
\Rightarrow\langle\vec{u}'\cdot\nabla\vec{T}'\rangle_h  \neq 0 ,
$$

and 

$$
\vec{\bar{u}}\cdot\nabla\vec{\bar{S}} = 0,\\
\Rightarrow\langle\vec{u}'\cdot\nabla\vec{\bar{S}}\rangle_h = \langle\vec{u}'\rangle_h\cdot\nabla\vec{\bar{S}} = 0,\\
\Rightarrow\langle\vec{\bar{u}}\cdot\nabla\vec{S}'\rangle_h = \vec{\bar{u}}\cdot\nabla\langle\vec{S}'\rangle_h = 0,\\
\Rightarrow\langle\vec{u}'\cdot\nabla\vec{S}'\rangle_h  \neq 0 ,
$$

The momentum equation can becomes

$$\partial_t (\vec{\bar{u}}+\vec{u}') + U_{bg}\partial_x(\vec{\bar{u}}+\vec{u}') + \partial_z U_{bg} (\bar{w}+w')\vec{e}_x + (\vec{\bar{u}}+\vec{u}')\cdot\nabla(\vec{\bar{u}}+\vec{u}')$$ 

$$= -\nabla p + \frac{Pr}{Pe}\nabla^2(\vec{\bar{u}}+\vec{u}') + \frac{4\pi^2 Ri}{R_\rho-1}((\bar{T}+T')-(\bar{S}+S'))\vec{e}_z,$$

$$\Leftrightarrow 
    \partial_t (\vec{\bar{u}}+\vec{u}')+ U_{bg}\partial_x\vec{u}' + \partial_z U_{bg}w'\vec{e}_{x}
    + \vec{\bar{u}}\cdot\nabla\vec{\bar{u}} + \vec{u}'\cdot\nabla\vec{\bar{u}} + \vec{\bar{u}}\cdot\nabla\vec{u}' + \vec{u}'\cdot\nabla\vec{u}'$$

$$= -\nabla p + \frac{Pr}{Pe}\nabla^2(\vec{\bar{u}}+\vec{u}') + \frac{4\pi^2 Ri}{R_\rho-1}((\bar{T}-\bar{S})+(T'-S'))\vec{e}_z,$$

then we horizontally average this equation, and get

$$
    \langle\partial_t (\vec{\bar{u}}+\vec{u}')\rangle_h
    + \langle U_{bg}\partial_x\vec{u}'\rangle_h + \langle\partial_z U_{bg}w'\vec{e}_x\rangle_h
    + \langle\vec{\bar{u}}\cdot\nabla\vec{\bar{u}}\rangle_h + \langle\vec{u}'\cdot\nabla\vec{\bar{u}}\rangle_h + \langle\vec{\bar{u}}\cdot\nabla\vec{u}'\rangle_h + \langle\vec{u}'\cdot\nabla\vec{u}'\rangle_h$$

$$= -\nabla \langle p\rangle_h + \langle\frac{Pr}{Pe}\nabla^2(\vec{\bar{u}}+\vec{u}')\rangle_h + \langle\frac{4\pi^2 Ri}{R_\rho-1}((\bar{T}-\bar{S})+(T'-S'))\vec{e}_z\rangle_h,$$

$$\Leftrightarrow 
    \partial_t \vec{\bar{u}}
    + \langle\vec{u}'\cdot\nabla\vec{u}'\rangle_h
    = -\nabla\langle p\rangle_h + \frac{Pr}{Pe}\nabla^2\vec{\bar{u}} + \frac{4\pi^2 Ri}{R_\rho-1}(\bar{T}-\bar{S})\vec{e}_z ,
$$

This equation can be called equation of the horizontally averaged velocity component
Now, we can take the equation of fluctuation = the original momentum equation - horizontally averaged equation

$$
    \partial_t \vec{u}'+ U_{bg}\partial_x\vec{u}' + \partial_z U_{bg}w'\vec{e}_x
    + \vec{u}'\cdot\nabla\vec{\bar{u}} + \vec{\bar{u}}\cdot\nabla\vec{u}' + \vec{u}'\cdot\nabla\vec{u}'
    - \langle\vec{u}'\cdot\nabla\vec{u}'\rangle_h
    = -\nabla p' + \frac{Pr}{Pe}\nabla^2\vec{u}' + \frac{4\pi^2 Ri}{R_\rho-1}(T'-S')\vec{e}_z,
$$

And,

$$\langle\vec{u}'\cdot\nabla\vec{u}'\rangle_h 
= \langle\begin{bmatrix}
    \vec{u}'\cdot\nabla u' \\
    \vec{u}'\cdot\nabla v' \\
    \vec{u}'\cdot\nabla w'
  \end{bmatrix}\rangle_h$$

$$= \langle\begin{bmatrix}
    u'\partial_x u'+v'\partial_y u'+w'\partial_z u'\\
    u'\partial_x v'+v'\partial_y v'+w'\partial_z v'\\
    u'\partial_x w'+v'\partial_y w'+w'\partial_z w'
  \end{bmatrix}\rangle_h$$

$$= \begin{bmatrix}
    \langle u'\partial_x u'+v'\partial_y u'+w'\partial_z u'\rangle_h \\
    \langle u'\partial_x v'+v'\partial_y v'+w'\partial_z v'\rangle_h \\
    \langle u'\partial_x w'+v'\partial_y w'+w'\partial_z w'\rangle_h 
  \end{bmatrix}$$

$$= \begin{bmatrix}
    \langle[\partial_x (u'u')-(\partial_x u')u']+[\partial_y (v'u')-\partial_y (v')u']+[\partial_z (w'u')-(\partial_z w')u']\rangle_h \\
    \langle[\partial_x (u'v')-(\partial_x u')v']+[\partial_y (v'v')-\partial_y (v')v']+[\partial_z (w'v')-(\partial_z w')v']\rangle_h \\
    \langle[\partial_x (u'w')-(\partial_x u')w']+[\partial_y (v'w')-\partial_y (v')w']+[\partial_z (w'w')-(\partial_z w')w']\rangle_h 
  \end{bmatrix}$$

$$= \begin{bmatrix}
    \langle[\partial_x (u'u')+\partial_y (v'u')+\partial_z (w'u')]-[\partial_x u'+\partial_y v'+\partial_z w']u'\rangle_h \\
    \langle[\partial_x (u'v')+\partial_y (v'v')+\partial_z (w'v')]-[\partial_x u'+\partial_y v'+\partial_z w']v'\rangle_h \\
    \langle[\partial_x (u'w')+\partial_y (v'w')+\partial_z (w'w')]-[\partial_x u'+\partial_y v'+\partial_z w']w'\rangle_h 
  \end{bmatrix}$$

$$= \begin{bmatrix}
    \langle\partial_x (u'u')+\partial_y (v'u')+\partial_z (w'u')\rangle_h \\
    \langle\partial_x (u'v')+\partial_y (v'v')+\partial_z (w'v')\rangle_h \\
    \langle\partial_x (u'w')+\partial_y (v'w')+\partial_z (w'w')\rangle_h 
  \end{bmatrix}
= \begin{bmatrix}
    \langle\partial_z (w'u')\rangle_h \\
    \langle\partial_z (w'v')\rangle_h \\
    \langle\partial_z (w'w')\rangle_h 
  \end{bmatrix}$$

$$\Leftrightarrow 
  \langle\vec{u}'\cdot\nabla\vec{u}'\rangle_h = \langle\partial_z (w'\vec{u}')\rangle_h
$$

Then, we have

$$
    \partial_t \vec{\bar{u}}
    + \partial_z \langle w'\vec{u}'\rangle_h
    = -\nabla\langle p\rangle_h + \frac{Pr}{Pe}\nabla^2\vec{\bar{u}} + \frac{4\pi^2 Ri}{R_\rho-1}(\bar{T}-\bar{S})\vec{e}_z,\\
    \partial_t \vec{u}'+ U_{bg}\partial_x\vec{u}' + \partial_z U_{bg}w'\vec{e}_x
    + \vec{u}'\cdot\nabla\vec{\bar{u}} + \vec{\bar{u}}\cdot\nabla\vec{u}' + \vec{u}'\cdot\nabla\vec{u}'
    - \partial_z \langle w'\vec{u}'\rangle_h
    = -\nabla (p-\langle p\rangle_h) + \frac{Pr}{Pe}\nabla^2\vec{u}' + \frac{4\pi^2 Ri}{R_\rho-1}(T'-S')\vec{e}_z,$$

The temperature equation:

$$\partial_t (\bar{T}+T') + U_{bg}\partial_x (\bar{T}+T') + (\vec{\bar{u}}+\vec{u}')\cdot\nabla (\bar{T}+T') - (\bar{w}+w') = \frac{1}{Pe}\nabla^2 (\bar{T}+T'),$$

$$\Leftrightarrow 
    \partial_t (\bar{T}+T') + U_{bg}\partial_x T' + \vec{\bar{u}}\cdot\nabla\bar{T} + \vec{u}'\cdot\nabla\bar{T} + \vec{\bar{u}}\cdot\nabla T' + \vec{u}'\cdot\nabla T'
    - (\bar{w}+w') = \frac{1}{Pe}\nabla^2 (\bar{T}+T'),$$

By using the same maner, we horizontal average above equation, and can get

$$
    \partial_t \bar{T} + \langle\vec{u}'\cdot\nabla T'\rangle_h - \bar{w} = \frac{1}{Pe}\nabla^2 \bar{T} = \frac{1}{Pe}\partial_z^2 \bar{T},
$$

then

$$
\partial_t T'  + U_{bg}\partial_x T' + \vec{u}'\cdot\nabla\bar{T} + \vec{\bar{u}}\cdot\nabla T' + \vec{u}'\cdot\nabla T' - \langle\vec{u}'\cdot\nabla T'\rangle_h - w' = \frac{1}{Pe}\nabla^2 T',
$$

due to $\vec{\bar{u}}\cdot\nabla\bar{T}=0$
And,

$$
\langle\vec{u}'\cdot\nabla T'\rangle_h = \partial_z \langle w'T' \rangle_h
$$

Then, we have

$$\partial_t \bar{T} + \partial_z \langle w'T' \rangle_h - \bar{w} = \frac{1}{Pe}\partial_z^2 \bar{T},$$

$$\partial_t T' + U_{bg}\partial_x T' + \vec{u}'\cdot\nabla\bar{T} + \vec{\bar{u}}\cdot\nabla T' + \vec{u}'\cdot\nabla T' - \partial_z \langle w'T' \rangle_h - w' = \frac{1}{Pe}\nabla^2 T',$$

The sanility equation:

$$\partial_t (\bar{S}+S') + U_{bg}\partial_x (\bar{S}+S') + (\vec{\bar{u}}+\vec{u}')\cdot\nabla (\bar{S}+S') - R_\rho (\bar{w}+w')= \frac{\tau}{Pe}\nabla^2 (\bar{S}+S'),$$

$$\Leftrightarrow 
    \partial_t (\bar{S}+S') + U_{bg}\partial_x S' + \vec{\bar{u}}\cdot\nabla\bar{S} + \vec{u}'\cdot\nabla\bar{S} + \vec{\bar{u}}\cdot\nabla S' + \vec{u}'\cdot\nabla S' - R_\rho w' = \frac{\tau}{Pe}\nabla^2 (\bar{S}+S'),$$

then

$$\partial_t \bar{S} + \langle\vec{u}'\cdot\nabla S'\rangle_h = \frac{\tau}{Pe}\nabla^2 \bar{S}=  \frac{\tau}{Pe}\partial_z^2 \bar{S},$$

$$\partial_t S' + \vec{u}'\cdot\nabla\bar{S} + \vec{\bar{u}}\cdot\nabla S' + \vec{u}'\cdot\nabla S' - \langle\vec{u}'\cdot\nabla S'\rangle_h - R_\rho w' = \frac{\tau}{Pe}\nabla^2 S',$$

due to $\vec{\bar{u}}\cdot\nabla\bar{S}=0$
And,

$$
\langle\vec{u}'\cdot\nabla S'\rangle_h = \partial_z \langle w'S' \rangle_h
$$

Then, we have

$$\partial_t \bar{S} + \partial_z \langle w'S'\rangle_h - R_\rho \bar{w}= \frac{\tau}{Pe}\partial_z^2 \bar{S},$$

$$\partial_t S' + U_{bg}\partial_x S' + \vec{u}'\cdot\nabla\bar{S} + \vec{\bar{u}}\cdot\nabla S' + \vec{u}'\cdot\nabla S' - \partial_z \langle w'S'\rangle_h - R_\rho w' = \frac{\tau}{Pe}\nabla^2 S',$$