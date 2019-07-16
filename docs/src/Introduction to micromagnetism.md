# <center>Computing energy barrier of magnetic system by micromagnetic method
<center>  吕伯尧 </center>
<center>SA18168036 </center>
## Introduction to micromagnetism
In ferromagnetic material, electrons within a small scale(~1nm) are strongly coupled by the exchange interaction, which causes a local magnetic moment. The opinion of micromagnetism is:  Assuming all magnetic moments are localized and arranged like a solid lattice,and the magnetic moments centered on the lattice points have a same fixed modular length and rotate freely under the influence of the effective magnetic field.
The magnetostatic energy and magnetic field has several terms of contributions:
##### Exchange energy

  ```math
  E_\mathrm{ex} = \int_{V} A (\nabla \vec{m})^2 \mathrm{d}V
  ```
  The exchange field is a sum of exchange interactions from 6 nearest neighbor cells:
```math
\vec{H}_{i, e x}=\frac{A}{\mu_{0}} \sum_{<i, j>} \vec{m}_{j}
```
##### Zeeman energy
  ```math
  E_\mathrm{ex} = \int_{V} \vec{H} \cdot \vec{m} \mathrm{d}V
  ```
  Applied zeeman field:
```math
\vec{H}_\mathrm{i,zeeman}=\vec{H}_\mathrm{ext}

```
##### Anisotropy

```math
E_\mathrm{anis} = -\int_{V} K_{u} (\vec{m} \cdot \hat{u})^2 \, dV
```
Set the angles between the anistropy axis and x,y,z coordinates as α,β,γ respectively, the anistropy field can be written as :
```math
\vec{H}_\mathrm{i,anis}=\frac{2 \mathcal{K}_{u}}{\mu_{s}}\left({m_{x}  \cos^{2}\alpha \cdot\vec{i}} +{m_{y}  \cos^{2}\beta \cdot \vec{j}}+ m_{z} \cos^{2}\gamma \cdot \vec{k} \right)

```
##### Bulk DMI energy

  ```math
  E_{\mathrm{dmi}} = \int_V D \vec{m} \cdot (\nabla \times \vec{m}) \, \mathrm{d}V
  ```


##### Demagnetization
```math
\mathscr{E}_{\text { demag }}=-\frac{1}{2} \overrightarrow{\mathbf{M}} \cdot \overrightarrow{\mathbf{B}}_{\text { demag }}
```
The demagnetization is from dipolar interactions between spins, the field is as the format:
```math
\overrightarrow{\mathbf{B}}_{\text { demag } }=\widehat{\mathbf{K}}_{i j} * \overrightarrow{\mathbf{M}}_{j}
```
In JuMag we use FFT to calculate the demag kernel $/hat{k_{ij}} $



  


  
