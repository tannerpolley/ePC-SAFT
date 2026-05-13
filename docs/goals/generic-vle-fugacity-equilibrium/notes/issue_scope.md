# Issue #90 Scope

- Issue URL: https://github.com/tannerpolley/ePC-SAFT/issues/90
- Title: Generic VLE/fugacity-equilibrium solver for volatile neutral species
- Current state: OPEN

## Scope

- direct volatile partial pressure from liquid fugacity where valid
- bubble pressure
- bubble temperature
- dew pressure
- dew temperature
- TP flash
- diagnostics report route used

## Out of Scope

- do not create MEA bubble-pressure-specific public APIs
- do not assume ions distribute to vapor unless the problem explicitly models it

## Policy Reminders

- keep `epcsaft` generic and application-neutral
- no finite difference
- use the approved derivative policy for explicit algebraic derivatives and solved states
