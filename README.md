


# SQIPrime

This is a proof-of-concept implementation of **SQIPrime2D**, as detailed in [SQIPrime](https://eprint.iacr.org/2024/773) by Duparc-Fouotsa, using [SageMath](https://www.sagemath.org).

This version is designed to demonstrate the validity of the mathematics behind SQIPrime.

Our code is partially based on the following works:

- [Theta-Sagemath](https://github.com/ThetaIsogenies/two-isogenies) by [Dartois-Maino-Pope-Robert](https://eprint.iacr.org/2023/1747).
- [QFESTA-SageMath](https://github.com/hiroshi-onuki/QFESTA-SageMath/tree/main) by [Hiroshi-Kohei](https://link.springer.com/chapter/10.1007/978-3-031-68388-6_4).
- [Learning-to-SQI](https://github.com/LearningToSQI/SQISign-SageMath) by [Corte-Real Santos-Pope](https://learningtosqi.github.io/).

**Note:** Significant modifications were required to efficiently implement SQIPrime. For instance, we adapted the code from [Theta-Sagemath](https://github.com/ThetaIsogenies/two-isogenies) to work over supersingular elliptic curves defined over $\mathbb{F}_{p^4}$. This allows for efficient evaluation of points  $P \in E$ whose order divides either $p-1$ or $p+1$.

## Usage

### Requirements

This implementation requires a recent version of SageMath (v10.3 and above), which can be downloaded [here](https://doc.sagemath.org/html/en/installation/index.html).

Additionally, we use the `SHAKE` function from the `pycryptodome` library, which can be installed using the following command in your terminal:

```
pip install -r pycryptodome
```


### Tests
To test SQIPrime, simply run 

```Terminal
sage SQIprime_test.sage
```
in your terminal. 

You can also perform a bechmark of the differnet fucntion of SQIPrime by running

```Terminal
sage benchmark.sage
```
in your terminal.


### Usage Example

SQIPrime can be used as follows:

```python
#import
from src.SQIPrime import SQIPrime

# Initialize parameters
p117 = 2**247 * 79 - 1
q117 = 168118140144706967996895604212334429
e117 = 247  
sqiprime = SQIPrime(p117, q117, e117)

# Key generation
sk, pk = sqiprime.KeyGen()

# Sign a message
message = bytes("SQIPrime is a simple signature mechanism", 'utf-8')
sign = sqiprime.Sign(sk, pk, message)

# Verification
ver = sqiprime.Verif(pk, message, sign)
assert ver
```











