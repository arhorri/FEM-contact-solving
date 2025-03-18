# FEM-contact-solving
This repository implements a 1D finite element analysis of an axially loaded bar with variable cross-section. It compares Lagrange multiplier and penalty methods to enforce a displacement constraint, using vectorized energy calculations for improved performance and providing visualization of force-displacement behavior.

# The Problem
Consider an elastic bar with Young’s modulus of 210 GPa and Poisson’s ratio of 0.3. The bar has a circular cross-section, with one end having a 20 mm diameter and the other a 40 mm diameter. Four equally spaced nodes are assumed along its 1‑meter length.

![image](https://github.com/user-attachments/assets/2b13c684-0637-442e-a2f5-54c309ad7009)

Write programs to solve the problem using both the Lagrange multiplier and penalty methods, assuming that a rigid wall is placed at a distance of 1/0 mm from the bar’s end, and that a force (less than the material’s yield force) is applied at the mid-node. Plot the force–displacement curve at the point of force application for both methods.
In the penalty method, plot the interference between the bar and the rigid wall as a function of the penalty stiffness coefficient and analyze the effect of this coefficient on the solution accuracy.

# Intruduction
The subject of contact mechanics focuses on understanding and predicting the interactions that occur when bodies come into contact. This presentation introduces the theoretical and numerical foundations of contact mechanics, beginning with a simple mass-spring system as a model to illustrate how contact forces arise when a mass engages with a rigid surface. Key topics include the formulation of energy variations and the application of different constraint methods, such as the Lagrange multiplier and penalty methods, to accurately simulate the behavior of contacting bodies.

In exploring these methods, the presentation highlights that the Lagrange multiplier technique ensures the strict fulfillment of the contact constraints by introducing additional variables, while the penalty method, though simpler, can result in slight penetrations between contacting surfaces. The discussion extends to finite element analysis (FEA), where the contact between two bars is modeled, underscoring the challenges that arise from the changing topology of the system matrix when contact constraints become active or inactive.

Overall, this introductory material lays the groundwork for more advanced studies in computational contact mechanics, emphasizing its crucial role in accurately modeling and analyzing interactions in engineering systems, particularly in contexts such as structural design and material behavior under load.

# Problem Statement and Geometry
### ###

We have a tapered (conical) bar:

Length: $L = 1$ m.

Diameter at one end (fully fixed end): $d_1 = 20$ mm

Diameter at the other end (near the rigid wall): $d_2 = 40$ mm

Material properties:

- Young’s modulus: $E = 210$ GPa
- Poisson’s ratio: $\nu = 0.3$ (not directly needed for a 1D axial analysis)

The bar is discretized into 4 finite elements (and thus 5 nodes). Nodes are equally spaced along the 1-meter length:

- Node 1: Fully fixed (clamped)
- Nodes 2, 3, 4: Interior nodes
- Node 5: Free end (initially with a gap of $g_0 = 0.1$ mm to a rigid wall)

A force $F$ (below the yield limit) is applied at Node 3 (the middle node).

#### Contact Modeling Between Node 5 and the Rigid Wall

If $u_5$ (displacement of Node 5) $< 0.1$ mm, no contact force.

If $u_5 \geq 0.1$ mm, contact is activated.

#### Solution Methods

The problem will be solved in two ways:

1. **Lagrange multiplier method** (strict constraint enforcement).
2. **Penalty method** (large spring stiffness $k_p$ to simulate contact).

We will compare the force–displacement curve at the load application point (Node 3) for both methods and examine how penalty stiffness $k_p$ affects the penetration $u_5 - g_0$.

### 1.1. Cross-Section Calculation

The diameter varies linearly from \$d\_1 = 20\$ mm to \$d\_2 = 40\$ mm over a length of 1 m. Converting to meters:

```math
d(x) = 0.02 + (0.04 - 0.02)x, \quad \text{where } x \in [0, 1] \text{ m}
```

At a position \$x\$:

```math
d(x) = 0.02 + 0.02x \text{ (meters)}
```

Cross-sectional area:

```math
A(x) = \pi \left( \frac{d(x)}{2} \right)^2 = \pi (0.01 + 0.01x)^2
```

The bar is divided into 4 equal elements, each of length 0.25 m. The cross-sectional area is evaluated at the midpoint of each element. For element \$e\$, the midpoint is:

```math
x_{\text{mid}, e} = \left(e - 0.5\right) \times \frac{L}{4}, \quad e = 1, 2, 3, 4
```

For example, for Element 1 (between Nodes 1 and 2), the midpoint is \$x = 0.125\$ m:

```math
d_{m,1} = 0.02 + 0.02 \times 0.125 = 0.0225 \text{ m}
```

```math
A_1 = \pi \times \left(\frac{0.0225}{2}\right)^2
```

Similarly, compute for elements 2, 3, and 4 using midpoints at \$x = 0.375\$, \$0.625\$, and \$0.875\$ m, respectively.

These concise formulas facilitate easy implementation and documentation in computational mechanics projects.



### 1.2. Elemental Strain
Calculate axial strain $\varepsilon_e$ for each element:

```math
\varepsilon_e = \frac{u_{e+1} - u_e}{L_e}
```

**Where:**
- $u_e$, $u_{e+1}$: Displacements at nodes $e$ and $e+1$
- $L_e$: Length of element $e$

### 1.3. Elemental Potential Energy
Compute potential energy $\Pi_e$ for each element:

```math
\Pi_e = \frac{1}{2} E A_e \varepsilon_e^2 L_e
```

**Where:**
- $E$: Young's modulus
- $A_e$: Cross-sectional area at element midpoint

### 1.4. Total Potential Energy
Sum elemental energies and subtract external work $F_3 u_3$:

```math
\Pi_{\text{total}} = \sum_{e=1}^{n} \Pi_e - F_3 u_3
```

**Where:**
- $F_3$: Force applied at Node 3
- $u_3$: Displacement at Node 3

These concise formulas facilitate easy implementation and documentation in computational mechanics projects.

### Lagrange Multiplier Method
#Here’s a revised version of your text incorporating the GitHub mathematical expression formatting:

---


# 2. Contact Constraint Methods: Lagrange Multiplier vs. Penalty Method

## 2.1. Lagrange Multiplier Method

### Extracted Formula and Conditions

The Lagrange multiplier method augments the original energy functional by adding a term that enforces the contact constraint. In its basic form, the augmented potential energy is written as:

```math
\Pi(u, \lambda) = \Pi_0(u) + \lambda \cdot (u - h)
```

where:

- $`\Pi_0(u)`$ is the original energy of the system (e.g., the energy of a mass-spring system),
- $`u`$ is the displacement,
- $`h`$ represents the gap or the prescribed contact position,
- $`\lambda`$ is the Lagrange multiplier, which is interpreted as the reaction force ($`R_N`$) when contact is active.

The method also requires satisfying the **Kuhn-Tucker conditions**:

- **Constraint Condition:** $`u - h \leq 0`$ (the displacement should not exceed the gap)
- **Non-Negativity:** $`\lambda \geq 0`$ (the reaction force cannot be negative)
- **Complementarity:** $`\lambda \cdot (u - h) = 0`$ (either there is no penetration, or there is a reaction force, but not both)

### Explanation

The Lagrange multiplier method enforces the contact condition exactly. When the constraint ($`u = h`$) is active, the Lagrange multiplier $`\lambda`$ becomes the reaction force that prevents penetration into the rigid support. This method introduces an additional unknown ($`\lambda`$) into the system, increasing the complexity of the equations but ensuring strict adherence to the contact constraint. As a result, it is highly accurate for problems requiring precise constraint enforcement.

---

## 2.2. Penalty Method

### Extracted Formula

The penalty method modifies the energy functional by adding a penalty term that penalizes any violation of the contact constraint. A typical form of the augmented energy functional in the penalty method is:

```math
\Pi(u) = \Pi_0(u) + \frac{1}{2} \cdot \varepsilon \cdot [\max(0, u - h)]^2
```

where:

- $`\varepsilon`$ is the penalty stiffness (a large number that controls the penalty severity for constraint violation),
- $`\max(0, u - h)`$ ensures the penalty applies only when the displacement $`u`$ exceeds $`h`$ (i.e., when penetration occurs).

The corresponding reaction force is approximated by:

```math
R_N = \varepsilon \cdot \max(0, u - h)
```

### Explanation

The penalty method enforces the contact constraint approximately. If the displacement $`u`$ exceeds $`h`$ (indicating penetration), the penalty term imposes a significant energy cost, pushing the solution back towards $`u = h`$. The reaction force $`R_N`$ is proportional to the penetration, scaled by the penalty stiffness $`\varepsilon`$. However, since the constraint is not enforced exactly, small penetration may persist if $`\varepsilon`$ is not sufficiently high.

The advantage of this approach is that it does not introduce additional unknowns into the system, preserving the size of the original equations. However, setting $`\varepsilon`$ too high can lead to numerical stiffness, while setting it too low may result in noticeable penetration. The choice of $`\varepsilon`$ is therefore crucial for balancing constraint enforcement and numerical stability.

---
