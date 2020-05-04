# **Hephaestus**
### Temperature Sensor Calculator

Makes available functions for:

* ITS-90 compatible temperature sensors (Standard Platinum Resistance Thermometers), allowing:
  * Retrieval of reference resistance ratio _W<sub>R</sub>_ from _t<sub>90R</sub>_ reference temperature
  * Retrieval of reference temperature _t<sub>90R</sub>_ from _W<sub>R</sub>_ reference resistance ratio
  * Given sensor coefficient's: _sub-range_, _R<sub>tpw</sub>_, _a_, _b_, _c_, _d_, _W<sub>660</sub>_, _c<sub>1</sub>_, _c<sub>2</sub>_, _c<sub>3</sub>_, _c<sub>4</sub>_, _c<sub>5</sub>_
    * Retrieval of  _t<sub>90</sub>_ temperature from sensor's resistance _res_ 
    * Retrieval of  _res_ resistance from sensor's _t<sub>90</sub>_ temperature. 
* IEC 60751 compatible temperature sensors (Platinum Resistance Thermometers), allowing:
  * Retrieval of reference resistance _res<sub>R</sub>_ from _t<sub>90R</sub>_ reference temperature
  * Retrieval of reference temperature _t<sub>90R</sub>_ from _res<sub>R</sub>_ reference resistance
  * Given sensor coefficient's: _R<sub>0</sub>_, _A_, _B_, _C_
    * Retrieval of sensor's resistance _res_ from _t<sub>90</sub>_ temperature
    * Retrieval of temperature _t<sub>90</sub>_ from sensor's resistance _res_
* IEC 60584 compatible temperature sensors (Thermocouples), allowing, for a given _type_:
  * Retrieval of reference voltage _emf<sub>R</sub>_ from _t<sub>90R</sub>_ reference temperature
  * Retrieval of reference temperature _t<sub>90R</sub>_ from reference voltage _emf<sub>R</sub>_
  * Given sensor coefficient's: _a_, _b_, _c_, _d_
    * Retrieval of sensor's voltage _emf_ from _t<sub>90</sub>_ temperature
    * Retrieval of temperature _t<sub>90</sub>_ from sensor's voltage _emf_


## Installation

Fill this with instructions specific to Python:


## Usage examples

**ITS-90**
```python
# Reference functions (ºC vs Ohm/Ohm)

# Corrected for specific sensor (ºC vs Ohm)
```

**IEC 60751**
```python
# Reference functions (ºC vs Ohm)

# Corrected for specific sensor
```

**IEC 60584**
```python
# Reference functions (ºC vs mV)

# Corrected for specific sensor
```


### Notes


## Contributing

1. Fork it ( https://github.com/[my-github-username]/hephaestus/fork )
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create a new Pull Request