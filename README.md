# Swirl Injector Design

---

### Coaxial_Swirl_Injector.m

- Input fluid parameters, mass flow rate, density, and pressure drop.
- Input spray cone angle of inner element, number of tangential inlets, and hydraulic loss coefficients.
- Input op_cl configurations for inner and outer elements and select coefficients of nozzle opening.
- Input inner element nozzle wall thickness.
- Input rec configurations to select mixing type, update mixing time for internal mixing selection.
- Run script and check outputs.
  
Common errors:
  
- Coefficient of nozzle opening is invalid for op_cl selections.
- Inner element is not accomodated within gas vortex of outer element.
- Spray angle of outer stage is too large for external mixing.

### Swirl_Injector.m

- Input fluid parameters, mass flow rate, density, and pressure drop.
- Input ranges or values for filling efficiency and coefficient of nozzle opening.
- Run script and check output figures.

### Swirl_Injector_Dynamic_Response.m

- Input injector dimensions, fluid density, and pressure drop.
- Select frequency range and artificial viscosity coefficient.
- Run script and check output figure.
