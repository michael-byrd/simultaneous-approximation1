# Author: Michael Byrd
# Date: 2024

from typing import Union, List, Dict, Tuple


class BivariatePolynomial:
    """
    A class to represent a bivariate polynomial.

    Attributes:
    -----------
    coefficients : dict
        A dictionary where keys are tuples representing the powers of x and y
        (in that order), and values are the corresponding coefficients.
    deg : int
        The degree of the polynomial.
    """

    def __init__(self, coefficients: Union[List[float], Dict[Tuple[int, int], float]]) -> None:
        """
        Initializes the BivariatePolynomial with either a list of coefficients
        or a dictionary of coefficients.

        Parameters:
        -----------
        coefficients : list or dict
            If a list is provided, it represents the coefficients of the polynomial
            in a specific order. The list should be in the format:
            [c00, c10, c01, c20, c11, c02, c30, c21, c12, c03, ...]
            where cij represents the coefficient of x^i * y^j.
            If a dictionary is provided, it should be in the format:
            {(i, j): cij, ...}
            where (i, j) is a tuple representing the powers of x and y, and cij is
            the coefficient of x^i * y^j.

        Example:
        --------
        For a polynomial c00 + c10*x + c01*y + c20*x^2 + c11*x*y + c02*y^2,
        the coefficients list should be [c00, c10, c01, c20, c11, c02].

        If a dictionary is provided, it could look like:
        {(0, 0): c00, (1, 0): c10, (0, 1): c01, (2, 0): c20, (1, 1): c11, (0, 2): c02}

        Notes:
        ------
        If the coefficients list is provided, the class will convert it into a
        dictionary format and initialize the polynomial. If a dictionary is provided,
        it directly initializes the polynomial.
        """
        self.deg: int = -1
        if type(coefficients) is list:
            # Initialize the coefficients dictionary
            self.coefficients: Dict[Tuple[int, int], float] = {}
            coefficient_list_length = len(coefficients)
            coefficient_list_index, degree_counter = 0, -1

            # Convert the coefficients list to a dictionary format
            while coefficient_list_index < coefficient_list_length:
                degree_counter += 1
                x_degree = degree_counter
                y_degree = 0

                # Iterate through the list and assign coefficients to the corresponding (x, y) powers
                while y_degree <= degree_counter and coefficient_list_index < coefficient_list_length:
                    if coefficients[coefficient_list_index] != 0:
                        self.coefficients[(x_degree, y_degree)] = coefficients[coefficient_list_index]
                    coefficient_list_index += 1
                    # For equal degree monomials (x^i * y^j), decrement x_degree and increment y_degree
                    x_degree -= 1
                    y_degree += 1
            self.calculate_degree()
        if type(coefficients) is dict:
            self.coefficients = coefficients
            self.reduce()

    def __str__(self) -> str:
        """
        Returns a string representation of the polynomial.
        :return: A string representation of the polynomial.
        """
        if not self.coefficients:
            return "0"

        terms = []

        for total_degree in range(self.deg + 1):
            for y_power in range(total_degree + 1):
                x_power = total_degree - y_power
                if (x_power, y_power) in self.coefficients:
                    coefficient = self.coefficients[(x_power, y_power)]

                    # Determine the sign and the absolute value of the coefficient
                    sign = "+" if coefficient > 0 else "-"
                    abs_coefficient = abs(coefficient)

                    # Format the coefficient: omit 1 or -1 unless the term is constant
                    if abs_coefficient == 1 and (x_power > 0 or y_power > 0):
                        formatted_coefficient = ""  # Omit 1 unless it's a constant term
                    else:
                        formatted_coefficient = str(abs_coefficient)

                    # Build the term
                    term = f"{formatted_coefficient}"
                    if x_power > 0:
                        term += f"x^{x_power}" if x_power > 1 else "x"
                    if y_power > 0:
                        term += f"y^{y_power}" if y_power > 1 else "y"

                    # Prepend the sign (handle first term separately)
                    if not terms:
                        term = term if coefficient > 0 else f"-{term}"
                    else:
                        term = f" {sign} {term}"

                    terms.append(term)

        return "".join(terms) or "0"

    def __repr__(self) -> str:
        return str(self)

    def _derivative(self, variable):
        """
        Returns the derivative of the polynomial with respect to the specified variable.
        :param variable: 0 for x, 1 for y
        :return: BivariatePolynomial object representing the derivative
        """
        derivative = {}
        for (x_power, y_power), coefficient in self.coefficients.items():
            if variable == 0 and x_power > 0:
                # Differentiate with respect to x
                new_powers = (x_power - 1, y_power)
                derivative[new_powers] = coefficient * x_power
            elif variable == 1 and y_power > 0:
                # Differentiate with respect to y
                new_powers = (x_power, y_power - 1)
                derivative[new_powers] = coefficient * y_power

        return BivariatePolynomial(derivative)

    def derivative(self, variable, derivative_order=1):
        """
        Returns the derivative of the polynomial with respect to the specified variable.
        :param variable: 0 for x, 1 for y
        :param derivative_order: The order of the derivative
        :return: BivariatePolynomial object representing the derivative
        """
        derivative_poly = self
        for _ in range(derivative_order):
            derivative_poly = derivative_poly._derivative(variable)
        return derivative_poly

    def __add__(self, other):
        """
        Adds two bivariate polynomials.
        :param other: The other polynomial to add
        :return: A new BivariatePolynomial object representing the sum of the two polynomials
        """
        # Create a new dictionary to store the result without altering the original polynomials
        polynomial_sum = self.coefficients.copy()

        # Add terms from the other polynomial
        for term, coefficient in other.coefficients.items():
            polynomial_sum[term] = polynomial_sum.get(term, 0) + coefficient

            # Optionally remove terms with zero coefficients (clean up zero terms)
            if polynomial_sum[term] == 0:
                del polynomial_sum[term]

        return BivariatePolynomial(polynomial_sum)

    __rad__ = __add__

    def __mul__(self, other):
        """
        Multiplies the polynomial by a scalar or another polynomial.
        :param other: The scalar or polynomial to multiply by
        :return: A new BivariatePolynomial object representing the product
        """
        polynomial_product = {}
        if isinstance(other, (int, float)):
            for term, coefficient in self.coefficients.items():
                polynomial_product[term] = coefficient * other
        elif isinstance(other, BivariatePolynomial):
            for (x1, y1), coefficient1 in self.coefficients.items():
                for (x2, y2), coefficient2 in other.coefficients.items():
                    new_term = (x1 + x2, y1 + y2)
                    polynomial_product[new_term] = polynomial_product.get(new_term, 0) + coefficient1 * coefficient2
        else:
            raise TypeError("Unsupported type for multiplication: {type(other)}")
        return BivariatePolynomial(polynomial_product)

    __rmul__ = __mul__

    def __sub__(self, other):
        """
        Subtracts another polynomial from this polynomial.
        :param other: The polynomial to subtract
        :return: A new BivariatePolynomial object representing the difference
        """
        return self + (-1 * other)

    def calculate_degree(self) -> None:
        """
        Calculates the degree of the polynomial based on the powers of x and y
        
        The (total) degree of a bivariate polynomial is the maximum sum of the
        powers of x and y in any term of the polynomial.

        This method is called internally to set the degree of the polynomial.
        """
        if self.coefficients:
            self.deg = max([sum(powers) for powers in self.coefficients.keys()])
        # If the coefficient dictionary is empty (i.e., the polynomial is zero), set the degree to -1
        else:
            self.deg = -1

    def remove_zero_coefficients(self) -> None:
        """
        Removes all items from the coefficients dictionary whose value is 0.
        """
        self.coefficients = {monomial: coefficient for monomial, coefficient in self.coefficients.items() if
                             coefficient != 0}

    def reduce(self) -> None:
        """
        Reduces the polynomial by removing zero coefficients and
        determining the polynomial's degree.

        This method is called internally if a dictionary of coefficients
        is provided during initialization.
        """
        self.remove_zero_coefficients()
        self.calculate_degree()

    def evaluate(self, point):
        """
        Evaluates the polynomial at a given point (x, y).
        :param point: A tuple (x, y) representing the point at which to evaluate the polynomial
        :return: The value of the polynomial at the given point
        """
        x, y = point
        value = 0
        for (x_power, y_power), coefficient in self.coefficients.items():
            value += coefficient * (x ** x_power) * (y ** y_power)
        return value

    def gradient(self):
        """
        Returns the gradient of the polynomial as a list of two BivariatePolynomial objects,
        :return: [dP/dx, dP/dy]
        """
        return [self.derivative(0), self.derivative(1)]
