###########
# Imports #
###########

from sympy.core.numbers import pi, E, oo
from math import floor, ceil
import math as mt
from scipy.special import polygamma, gammainc, gammaincc, erf, erfc
import scipy.special as ss
from rich import print
from rich.progress import (
    Progress,
    SpinnerColumn,
    BarColumn,
    MofNCompleteColumn,
    TimeElapsedColumn,
    TextColumn,
)
from prompt_toolkit import prompt
from prompt_toolkit.key_binding import KeyBindings
from numpy import (
    linspace,
    exp,
    log,
    sqrt,
    abs,
    sin,
    cos,
    tan,
    arcsin,
    arccos,
    arctan,
    sinh,
    cosh,
    tanh,
    arcsinh,
    arccosh,
    arctanh,
    array,
    cross,
)
from numpy.linalg import norm, det

# Sympy uses symbolic Pi and e, which cannot be graphed by matplotlib
# so later np.pi np.e are explicitly used instead.
import numpy as np
from sympy import (
    diff,
    idiff,
    Integral,
    integrate,
    limit,
    solve,
    linsolve,
    pprint,
    simplify,
    symbols,
    gamma,
    lowergamma,
    uppergamma,
)
import matplotlib.pyplot as plt
import datetime, logging, random, readline, os, time, warnings


###########
# Presets #
###########

warnings.filterwarnings("ignore")

readline.parse_and_bind("tab: complete")

x, y, z = symbols("x, y, z")
history = []  # Stores calculation history

input = prompt

key_bindings = KeyBindings()

idx = len(history)


# Custom up/down key bindings


@key_bindings.add("up")
def _(event):
    global idx

    if idx != -1:
        idx -= 1
        buffer = event.app.current_buffer
        buffer.text = history[idx]
        buffer.cursor_position = len(buffer.text)
        idx = idx % len(history)
    else:
        idx = 0
        buffer = event.app.current_buffer
        buffer.text = history[idx]
        buffer.cursor_position = len(buffer.text)


@key_bindings.add("down")
def _(event):
    global idx

    if idx != len(history):
        idx += 1
        idx = idx % len(history)
        buffer = event.app.current_buffer
        buffer.text = history[idx]
        buffer.cursor_position = len(buffer.text)

    else:
        idx = len(history) - 1
        buffer = event.app.current_buffer
        buffer.text = history[idx]
        buffer.cursor_position = len(buffer.text)


#####################
# Calculation funcs #
#####################


def simp():
    print(
        '\n[bold bright_green](Current Screen: Simple Calculation Screen)[/bold bright_green]\n("q" to quit)\nEnter any expression to start:\n'
    )

    while True:
        expr = input(key_bindings=key_bindings)

        try:
            if expr != "q":
                # Easy to exit unlike VIM ;)
                result = eval(replace_graph(expr))
                history.append(str(result))
                print("\n[bright_yellow]Result: [/bright_yellow]", end="")
                pprint(result)
                print()

            else:
                break

        except:
            print()
            logging.error(f'Could not parse: "{expr}"\n')


def derive(function: str, order: str):
    calc = replace_expr(function)

    if check_order(order) is True:
        global df
        df = diff(calc, x, order)

        print(
            f"\nDerivative of [bright_magenta]{replace_graph(function)}[/bright_magenta] with order {order} is:\n"
        )
        print_expr(df)

        df = replace_expr(str(df))
        check_simp(df)
        history.append(df)

        df = replace_graph(df)

        print("\n[bright_yellow]Show graph of area? (y/n)[/bright_yellow]")
        show = input("(Exit the graph window when you are finished to continue) ")

        if show == "y":
            try:
                print()
                with Progress(
                    SpinnerColumn(finished_text="[bright_green]√[/bright_green]"),
                    TextColumn("[bright_yellow]Loading graph...[/bright_yellow]"),
                    BarColumn(),
                    TimeElapsedColumn(),
                    transient=False,
                ) as progress:
                    task = progress.add_task("", total=100)

                    while not progress.finished:
                        progress.update(task, advance=2)
                        time.sleep(random.randint(2, 5) / 1000)

                x_arr = linspace(-50, 50, 200000)

                def f(x: np.array):
                    return eval(replace_graph(calc))

                def dif(x: np.array):
                    return eval(df)

                title = "Function (red) and derivative (blue)"
                plt.title(title)
                plt.xlabel("x", weight="bold")
                plt.ylabel("y", rotation=0, weight="bold")
                plt.plot(x_arr, f(x_arr), color="red", label="Function")
                plt.plot(x_arr, dif(x_arr), color="blue", label="Derivative")
                plt.axis([-7.5, 7.5, -7.5, 7.5])
                plt.legend(loc="lower left")
                plt.grid()
                plt.show()

                print("\nExited graph.\n")

            except:
                plt.close()
                print("\n")
                logging.warning("Could not graph function.")
                print("\nExited graph.")


def partial_derive(function: str, var: str, order: str):
    calc = replace_expr(function)

    if check_order(order) is True:
        df = diff(calc, var, order)

        print(
            f"\nPartial derivative of [bright_magenta]{replace_graph(function)}[/bright_magenta] in respect to {var} of order {order} is:\n"
        )
        print_expr(df)

        df = replace_expr(str(df))
        check_simp(df)
        history.append(df)


def implicit_derive(circ: str, order: str):
    calc = replace_expr(circ)
    left = eval(calc[: calc.find("=")])
    right = eval(calc[calc.find("=") + 1 :])

    if check_order(order) is True:
        df = idiff(left - right, y, x, order)

        print(
            f"\nDerivative of [bright_magenta]{circ}[/bright_magenta] with order {order} is:\n"
        )
        print_expr(df)

        df = replace_expr(str(df))
        check_simp(df)
        history.append(df)


def antiderive(function: str):
    calc = replace_expr(function)
    F = Integral(calc, x).doit()

    if "Integral" in str(F):
        logging.warning("Cannot compute integral.\n")
        return ""

    print(
        f"\nAntiderivative of [bright_magenta]{replace_graph(function)}[/bright_magenta] is:\n"
    )
    print_expr(F)

    F = replace_expr(str(F))
    check_simp(F)
    history.append(F)

    F = replace_graph(F)

    print("\n[bold]Don't forget to add a constant![/bold]\n")

    print("\n[bright_yellow]Show graph of area? (y/n)[/bright_yellow]")
    show = input("(Exit the graph window when you are finished to continue) ")

    if show == "y":
        try:
            print()
            with Progress(
                SpinnerColumn(finished_text="[bright_green]√[/bright_green]"),
                TextColumn("[bright_yellow]Loading graph...[/bright_yellow]"),
                BarColumn(),
                TimeElapsedColumn(),
                transient=False,
            ) as progress:
                task = progress.add_task("", total=100)

                while not progress.finished:
                    progress.update(task, advance=2)
                    time.sleep(random.randint(2, 5) / 1000)

            x_arr = linspace(-100, 100, 200000)

            def f(x: np.array):
                return eval(replace_graph(calc))

            def af(x: np.array):
                return eval(F)

            title = "Function (red) and antiderivative (blue, C = 0)"
            plt.title(title)
            plt.xlabel("x", weight="bold")
            plt.ylabel("y", rotation=0, weight="bold")
            plt.plot(x_arr, f(x_arr), color="red", label="Function")
            plt.plot(x_arr, af(x_arr), color="blue", label="Antiderivative")
            plt.axis([-7.5, 7.5, -7.5, 7.5])
            plt.legend(loc="lower left")

            plt.grid()
            plt.show()

            print("\nExited graph.\n")

        except:
            plt.close()
            print("\n")
            logging.warning("Could not graph function.")
            print("\nExited graph.\n")


def def_int(function: str, low: str, up: str):
    calc = replace_expr(function)

    check_bound(low)

    calc_low = eval(replace_expr(low))

    check_bound(up)

    calc_up = eval(replace_expr(up))

    result = integrate(calc, (x, calc_low, calc_up))

    graph_up = eval(replace_expr(str(calc_up)))
    graph_low = eval(replace_expr(str(calc_low)))

    num_result = integrate(calc, (x, graph_low, graph_up)).evalf()
    # Composite functions usually do not have primitive antiderivatives
    # so calc-ultra is equipped with both symbolic and numerical answers.

    if (
        (str(result) == "nan")
        or ("I" in str(result))
        and ("Integral" not in str(result))
    ):
        logging.warning("Cannot compute integral because integral does not converge.")
        return ""

    if "Integral" not in str(result):
        print(
            f"\nCalculated integral of [bright_magenta]{replace_graph(function)}[/bright_magenta] from {low} to {up}. Final area is:\n"
        )
        print_expr(result)
        result = replace_expr(str(result))
        history.append(result)

    else:
        print("\nCannot express result symbolically.")

    nprint("\nNumeral evaluation/approximation:\n")
    print_expr(num_result)

    print("\n[bright_yellow]Show graph of area? (y/n)[/bright_yellow]")

    show = input("(Exit the graph window when you are finished to continue) ")

    if show == "y":
        try:
            print()
            with Progress(
                SpinnerColumn(finished_text="[bright_green]√[/bright_green]"),
                TextColumn("[bright_yellow]Loading graph...[/bright_yellow]"),
                BarColumn(),
                TimeElapsedColumn(),
                transient=False,
            ) as progress:
                task = progress.add_task("", total=100)

                while not progress.finished:
                    progress.update(task, advance=2)
                    time.sleep(random.randint(2, 5) / 1000)

            x_arr = linspace(-100, 100, 200000)

            if "log" in function:
                x_arr = linspace(
                    0.0001,
                    100,
                    200000,
                )

            # The factorial is just the gamma function shifted one unit left
            # Of course, the factorial is only defined for x >= 0

            if "factorial" in function:
                x_arr = linspace(
                    0,
                    100,
                    200000,
                )

            def f(x: np.array):
                return eval(replace_graph(calc))

            title = "Shaded area beneath function"
            plt.title(title)
            plt.xlabel("x", weight="bold")
            plt.ylabel("y", rotation=0, weight="bold")
            plt.plot(x_arr, f(x_arr), color="red", label="Function")

            plt.fill_between(
                x_arr,
                f(x_arr),
                where=[(x_arr > calc_low) and (x_arr < calc_up) for x_arr in x_arr],
                color="blue",
            )

            try:
                if graph_option == "f":
                    plt.axis([-7.5, 7.5, -7.5, 7.5])

                elif graph_option == "a":
                    # Adjusted graph view is sometimes better for
                    # large graphs with large bounds.
                    if (float(f(graph_low)) != 0) and (float(f(graph_up)) != 0):
                        plt.axis(
                            [
                                graph_low - 5,
                                graph_up + 5,
                                float(f(round(graph_low)))
                                - (
                                    float(f(round(graph_low)))
                                    + float(f(round(graph_up)))
                                )
                                / 2
                                - 1,
                                float(f(round(graph_up)))
                                + (
                                    float(f(round(graph_low)))
                                    + float(f(round(graph_up)))
                                )
                                / 2
                                + 1,
                            ]
                        )

                    elif (float(f(graph_low)) == 0) or (float(f(graph_up)) == 0):
                        plt.axis(
                            [
                                graph_low - 5,
                                graph_up + 5,
                                -(graph_up - graph_low) / 2,
                                (graph_up + graph_low) / 2,
                            ]
                        )

            except:
                plt.axis([-7.5, 7.5, -7.5, 7.5])

            plt.legend(loc="lower left")
            plt.grid()
            plt.show()

            return "\nExited graph.\n"

        except:
            plt.close()
            print("\n")
            logging.warning("Could not graph function.")
            return "\nExited graph.\n"

    else:
        return "\n[bright_yellow]Exiting Definite Integral Screen ... ... ...[/bright_yellow]\n"


def improp_int(function: str, low: str, up: str):
    function = replace_expr(function)

    check_bound(low)

    clow = eval(replace_expr(low))

    check_bound(up)

    cup = eval(replace_expr(up))

    try:
        improper_area = Integral(
            function, (x, clow, cup)
        ).principal_value()  # Cauchy Principal Value

        if "Integral" not in str(improper_area):
            print(
                f"\nCalculated improper integral of [bright_magenta]{function}[/bright_magenta] from {low} to {up}. Final area is:\n"
            )
            print_expr(improper_area)
            improper_area = replace_expr(str(improper_area))
            history.append(improper_area)

        else:
            print("Cannot compute improper integral.")

    except ValueError:
        improper_area = integrate(function, (x, clow, cup))

        if "Integral" not in str(improper_area):
            print(
                f"\nCalculated improper integral of [bright_magenta]{function}[/bright_magenta] from {low} to {up}. Final area is:\n"
            )
            print_expr(improper_area)
            improper_area = replace_expr(str(improper_area))
            history.append(improper_area)

        else:
            print("Cannot compute improper integral.")

    print()


def double_int(function: str, out_low: str, out_up: str, in_low: str, in_up: str):
    function = replace_expr(function)

    check_bound(out_low)

    out_low = eval(replace_expr(out_low))

    check_bound(out_up)

    out_up = eval(replace_expr(out_up))

    in_low = eval(replace_expr(in_low))

    in_up = eval(replace_expr(in_up))

    out_up = eval(replace_expr(str(out_up)))
    out_low = eval(replace_expr(str(out_low)))
    in_up = eval(replace_expr(str(in_up)))
    in_low = eval(replace_expr(str(in_low)))

    result = integrate(function, (y, in_low, in_up), (x, out_low, out_up))

    if "Integral" in str(result):
        logging.warning("Cannot compute integral.\n")
        return ""

    print(
        f"\nDouble integral of [bright_magenta]{function}[/bright_magenta] with inner bounds [bright_cyan]{in_low}[/bright_cyan] and [bright_cyan]{in_up}[/bright_cyan] and outer bounds {out_low} and {out_up} is:\n"
    )
    print_expr(result)
    result = replace_expr(str(result))
    history.append(result)
    print("")


def lim(expr: str, value: str):
    expr = replace_expr(expr)

    check_bound(value)

    value = float(eval(replace_expr(value)))

    l = limit(expr, x, value)

    if "Limit" in str(l):
        logging.warning("Cannot compute limit.")
        return ""

    if limit(expr, x, value, "+") != limit(expr, x, value, "-"):
        print(
            "\nThe limit does not exist (i.e., the limit approaching from the right does not equal the limit approaching from the left)."
        )
        print(
            "Use the one-side limit calculation of LimCalc to calculate the limit at one side.\n"
        )

    else:
        print(
            f"\nLimit of [bright_magenta]{expr}[/bright_magenta] as x approaches {value} is:\n"
        )
        print_expr(l)
        l = replace_expr(str(l))
        history.append(l)


def side_lim(expr: str, value: str, direction: str):
    expr = replace_expr(expr)

    check_bound(value)

    value = float(eval(replace_expr(value)))

    if direction == "left" or direction == "Left":
        direction_sign = "-"

    elif direction == "right" or direction == "Right":
        direction_sign = "+"

    else:
        print()
        logging.error("\nTypeError: Direction is neither right or left.")
        return ""

    l = limit(expr, x, value, dir=direction_sign)

    if "Limit" in str(l):
        logging.warning("\nCannot compute limit.")
        return ""

    print(
        f"\nLimit of [bright_magenta]{expr}[/bright_magenta] as x approaches {value} from the [bright_cyan]{direction}[/bright_cyan] is:\n"
    )
    print_expr(l)
    l = replace_expr(str(l))
    history.append(l)


def eq_solve(mode: int):
    eq_list = []

    if mode == "1":
        eq1 = input("\nEnter equation: ")
        left = eval(eq1[: eq1.find("=")])
        right = eval(eq1[eq1.find("=") + 1 :])
        eq_set = solve(left - right)
        if len(eq_set) == 0:
            print()
            logging.error("UnknownError: Cannot solve equation")
        else:
            print("\nx:\n")
            for i in range(0, len(eq_set)):
                pprint(eq_set[i])
                history.append(replace_expr(str(eq_set[i])))
                print()

    elif mode == "2":
        eq1 = input("Enter first equation: ")
        left1 = eval(eq1[: eq1.find("=")])
        right1 = eval(eq1[eq1.find("=") + 1 :])

        eq2 = input("Enter second equation: ")
        left2 = eval(eq2[: eq2.find("=")])
        right2 = eval(eq2[eq2.find("=") + 1 :])

        eq_set = str(linsolve((left1 - right1, left2 - right2), (x, y), (-1, 1)))
        eqs = eq_set.strip("{").strip("}").strip("(").strip(")")
        for value in eqs:
            eq_list.append(value)

        result = "".join(eq_list)
        print("\nx, y:\n")
        pprint(result)

    elif mode == "3":
        eq1 = input("Enter equation 1: ")
        left1 = eval(eq1[: eq1.find("=")])
        right1 = eval(eq1[eq1.find("=") + 1 :])

        eq2 = input("Enter equation 2: ")
        left2 = eval(eq2[: eq2.find("=")])
        right2 = eval(eq2[eq2.find("=") + 1 :])

        eq3 = input("Enter equation 3: ")
        left3 = eval(eq3[: eq3.find("=")])
        right3 = eval(eq3[eq3.find("=") + 1 :])

        eq_set = str(
            linsolve(
                (left1 - right1, left2 - right2, left3 - right3), (x, y, z), (-1, 1)
            )
        )
        eqs = eq_set.strip("{").strip("}").strip("(").strip(")")
        for value in eqs:
            eq_list.append(value)

        result = "".join(eq_list)
        print("\nx, y, z:\n")
        pprint(result)


################
# Helper funcs #
################


# These functions sanatize inputs and prevent errors


def check_simp(expr: str) -> bool:
    if str(simplify(eval(expr), evaluate=False)) != expr:
        print("\nSimplify/rewrite:\n")
        print_expr(simplify(expr, evaluate=False))

    else:
        return False


def check_order(order: str) -> bool:
    if ("." in order) or (order.isnumeric() == False):
        print()
        logging.error(
            "OrderError: Order of derivative calculation is not a valid number."
        )
        return False

    else:
        return True


def check_bound(bound: str):
    if (
        (bound.isnumeric() is False)
        and ("pi" not in bound)
        and ("e" not in bound)
        and ("E" not in bound)
        and ("-" not in bound)
        and ("." not in bound)
        and ("sqrt" not in bound)
        and ("oo" not in bound)
        and ("/" not in bound)
    ):
        print()
        logging.error("TypeError: Integration bound is not a number.")


def replace_expr(expr: str) -> str:
    expr = expr.strip(" ")

    if "^" in expr:
        expr = expr.replace("^", "**")

    # The extra 1 at the end is to recognize whether the "e"
    # is part of an expression or a constant on it's own.

    if "E" in expr:
        expr = expr.replace("E", "E1")

    if "e" in expr:
        expr = expr.replace("e", "E1")

    if "E1**" in expr:
        expr = expr.replace("E1**", "exp")

    if "E1xp" in expr:
        expr = expr.replace("E1xp", "exp")

    if "E1rf" in expr:
        expr = expr.replace("E1rf", "erf")

    if "pE1r" in expr:
        expr = expr.replace("pE1r", "per")

    if "wE1r" in expr:
        expr = expr.replace("wE1r", "wer")

    if "ln" in expr:
        expr = expr.replace("ln", "log")

    if "arc" in expr:
        expr = expr.replace("arc", "a")

    if "abs" in expr:
        expr = expr.replace("abs", "Abs")

    if "pi" in expr:
        expr = expr.replace("pi", str(np.pi))

    if "E1" in expr:
        expr = expr.replace("E1", str(np.e))

    return expr


def factorial(x):
    # Define our own factorial
    return ss.gamma(x + 1)


def replace_graph(function: str) -> str:
    # Sympy and Numpy trig functions uses different prefixes.
    # Thus this replacement algorithm is needed for graphing.
    if "asin" in function:
        function = function.replace("asin", "arcsin")

    if "acos" in function:
        function = function.replace("acos", "arccos")

    if "atan" in function:
        function = function.replace("atan", "arctan")

    if "asinh" in function:
        function = function.replace("asinh", "arcsinh")

    if "acosh" in function:
        function = function.replace("acosh", "arccosh")

    if "atanh" in function:
        function = function.replace("atanh", "arctanh")

    # Implement our own csc, sec, cot, csch, sech, and coth

    if "csc" in function:
        function = function.replace("csc", "1/sin")

    if "sec" in function:
        function = function.replace("sec", "1/cos")

    if "cot" in function:
        function = function.replace("cot", "1/tan")

    if "csch" in function:
        function = function.replace("csch", "1/sinh")

    if "sech" in function:
        function = function.replace("sech", "1/cosh")

    if "coth" in function:
        function = function.replace("coth", "1/tanh")

    if "Abs" in function:
        function = function.replace("Abs", "fabs")

    if "log" in function and ",x)" in function:
        function = function.replace("log", "mt.log")

    # Graph constant functions

    if "x" not in function:
        function = "0 * x + " + function

    if "gamma" in function:
        function = function.replace("gamma", "ss.gamma")

    if "polyss.gamma" in function:
        function = function.replace("polyss.gamma", "polygamma")

    if "lowergamma" in function:
        function = function.replace("lowergamma", "gammainc")

    if "uppergamma" in function:
        function = function.replace("uppergamma", "gammaincc")

    return function


def print_expr(text: str):
    printing_methods = {"p": lambda t: pprint(text), "n": lambda t: print(text)}

    try:
        printing_methods[print_option](text)

    except NameError:
        printing_methods["p"](text)


def nprint(text: str):
    print(text)
    time.sleep(0.04)


####################
# Driving programs #
####################


def main():
    with Progress(
        SpinnerColumn(finished_text="[bright_green]√[/bright_green]"),
        TextColumn("[green]Handling imports[/green]..."),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        transient=False,
    ) as progress:
        task = progress.add_task("", total=50)

        while not progress.finished:
            progress.update(task, advance=1)
            time.sleep(random.randint(2, 5) / 100)

    global x, y, z, history

    while True:
        instruct_path = (
            os.path.dirname(os.path.abspath(__file__))
            + "/texts/main_screen.txt"
            # TODO: make the PATH compatible with Windows
        )
        main = open(instruct_path, mode="r")
        for line in main.readlines():
            line = line.rstrip()
            if ("---" in line) or ("|" in line):
                nprint("[gold1]" + line + " [/gold1]")
            else:
                nprint(line)

        try:
            if date_option == "1":
                now = (datetime.datetime.now()).strftime("%Y/%m/%d")
            elif date_option == "2":
                now = (datetime.datetime.now()).strftime("%Y/%m/%d %H:%M:%S")
        except:
            now = (datetime.datetime.now()).strftime("%Y/%m/%d %H:%M:%S")

        print(f"\n(Time now is: {now})\n")
        print("[bold bright_green](Current Screen: Main Screen)[/bold bright_green]\n")
        print("[bright_purple]Enter Command: [/bright_purple]", end="")
        cmd = input()

        if cmd == "1":
            simp()

        elif cmd == "2":
            derivacalc()

        elif cmd == "3":
            intecalc()

        elif cmd == "4":
            limcalc()

        elif cmd == "5":
            algcalc()

        elif cmd == "6":
            settings()

        elif cmd == "7":
            nprint("\n[bright_yellow]Exiting Calc-ULTRA ... ... ...[/bright_yellow]\n")
            break

        else:
            print("\n")
            logging.warning(f'Invalid command:"{cmd}"\n')


"""
If you find this message, type 'hi' in the general discussions - sudoer-Huatao
"""


def derivacalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/derivacalc.txt"
    derivacalc = open(instruct_path, mode="r")
    for line in derivacalc.readlines():
        line = line.rstrip()
        if ("---" in line) or ("|" in line):
            nprint("[gold1]" + line + " [/gold1]")
        else:
            nprint(line)
    print()

    while True:
        try:
            nprint(
                "\n[bold bright_green](Current Screen: DerivaCalc Main Screen)[/bold bright_green]\n"
            )
            print("[bright_purple]Enter Command: [/bright_purple]", end="")
            cmd = input()

            if cmd == "1":
                nprint(
                    "\n[bold bright_green](Current Screen: Derivative Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ", key_bindings=key_bindings)
                order = input(
                    "Enter order of derivative calculation: ", key_bindings=key_bindings
                )
                derive(function, order)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current Screen: Partial Derivative Screen)[/bold bright_green]\n"
                )
                function = input(
                    "Enter a function containing x and y or x and y and z: ",
                    key_bindings=key_bindings,
                )
                var = input(
                    "Enter variable to differentiate in respect to: ",
                    key_bindings=key_bindings,
                )
                if var != "x" and var != "y" and var != "z":
                    print()
                    logging.error("Variable to differentite in respect to is invalid.")
                else:
                    order = input(
                        "Enter the order of partial derivative calculation: ",
                        key_bindings=key_bindings,
                    )
                    partial_derive(function, var, order)

            elif cmd == "3":
                nprint(
                    "\n[bold bright_green](Current Screen: Implicit Derivative Screen)[/bold bright_green]\n"
                )
                circ = input(
                    "Enter an equation containing x and y:", key_bindings=key_bindings
                )
                order = input(
                    "Enter order of implicit derivative calculation: ",
                    key_bindings=key_bindings,
                )
                implicit_derive(circ, order)

            elif cmd == "4":
                print("\n[bright_yellow]Exiting DerivaCalc ... ... ...[/bright_yellow]")
                break

            else:
                print("\n")
                logging.warning(f'Invalid command:"{cmd}"')
        except:
            print("\n")
            logging.error("UnknownError: An unknown error occured.")
            nprint("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def intecalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/intecalc.txt"
    intecalc = open(instruct_path, mode="r")
    for line in intecalc.readlines():
        line = line.rstrip()
        if ("---" in line) or ("|" in line):
            nprint("[gold1]" + line + " [/gold1]")
        else:
            nprint(line)
    print()

    while True:
        try:
            nprint(
                "[bold bright_green](Current Screen: InteCalc Main Screen)[/bold bright_green]\n"
            )
            print("[bright_purple]Enter Command: [/bright_purple]", end="")
            cmd = input()

            if cmd == "1":
                nprint(
                    "\n[bold bright_green](Current Screen: Antiderivative Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ", key_bindings=key_bindings)
                antiderive(function)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current Screen: Definite Integral Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ", key_bindings=key_bindings)
                lower_bound = input(
                    "\nEnter the lower bound: ", key_bindings=key_bindings
                )
                upper_bound = input(
                    "Enter the upper bound: ", key_bindings=key_bindings
                )
                print(def_int(function, lower_bound, upper_bound))

            elif cmd == "3":
                nprint(
                    "\n[bold bright_green](Current Screen: Improper Integral Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ", key_bindings=key_bindings)
                lower_bound = input(
                    "\nEnter the lower bound: ", key_bindings=key_bindings
                )
                upper_bound = input(
                    "Enter the upper bound: ", key_bindings=key_bindings
                )
                improp_int(function, lower_bound, upper_bound)

            elif cmd == "4":
                nprint(
                    "\n[bold bright_green](Current Screen: Double Integral Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ")
                outer_low = input(
                    "\nEnter the lower outer bound: ", key_bindings=key_bindings
                )
                outer_up = input(
                    "Enter the upper outer bound: ", key_bindings=key_bindings
                )
                inner_low = input(
                    "\nEnter the lower inner bound: ", key_bindings=key_bindings
                )
                inner_up = input(
                    "Enter the upper inner bound: ", key_bindings=key_bindings
                )
                double_int(function, outer_low, outer_up, inner_low, inner_up)

            elif cmd == "5":
                print("\n[bright_yellow]Exiting InteCalc ... ... ...[/bright_yellow]")
                break

            else:
                print("\n")
                logging.warning(f'Invalid command: "{cmd}"')
        except:
            print("\n")
            logging.error("UnknownError: An unknown error occured.")
            nprint("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def limcalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/limcalc.txt"
    limcalc = open(instruct_path, mode="r")
    for line in limcalc.readlines():
        line = line.rstrip()
        if ("---" in line) or ("|" in line):
            nprint("[gold1]" + line + " [/gold1]")
        else:
            nprint(line)
    print()

    while True:
        try:
            nprint(
                "\n[bold bright_green](Current Screen: LimCalc Main Screen)[/bold bright_green]\n"
            )
            print("[bright_purple]Enter Command: [/bright_purple]", end="")
            cmd = input()

            if cmd == "1":
                nprint(
                    "\n[bold bright_green](Current screen: Limit Screen)[/bold bright_green]\n"
                )
                expr = input("Enter an expression: ", key_bindings=key_bindings)
                value = input("Enter point of evaluation: ", key_bindings=key_bindings)
                lim(expr, value)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current screen: One-sided Limit Screen)[/bold bright_green]\n"
                )
                expr = input("Enter an expression: ", key_bindings=key_bindings)
                value = input("Enter point of evaluation: ", key_bindings=key_bindings)
                direction = input("Enter direction of limit ('left' or 'right'): ")
                side_lim(expr, value, direction)

            elif cmd == "3":
                print("\n[bright_yellow]Exiting LimCalc ... ... ...[/bright_yellow]")
                break

            else:
                print("\n")
                logging.warning(f'Invalid command: "{cmd}"')

        except:
            print("\n")
            logging.error("UnknownError: An unknown error occured.")
            nprint("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def algcalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/algcalc.txt"
    algcalc = open(instruct_path, mode="r")
    for line in algcalc.readlines():
        line = line.rstrip()
        if ("---" in line) or ("|" in line):
            nprint("[gold1]" + line + " [/gold1]")
        else:
            nprint(line)
    print()

    while True:
        try:
            nprint(
                "\n[bold bright_green](Current Screen: AlgCalc Main Screen)[/bold bright_green]\n"
            )
            print("[bright_purple]Enter Command: [/bright_purple]", end="")
            cmd = input()

            if cmd == "1":
                nprint(
                    "\n[bold bright_green](Current screen: Equation Solver Screen)[/bold bright_green]\n"
                )
                mode = ""

                while mode != "q":
                    print(
                        'Enter mode: 1 for one set equation, 2 for two, and 3 for three ("q" to quit): ',
                        end="",
                    )
                    mode = input()
                    eq_solve(mode)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current screen: Vector Calculation Screen)[/bold bright_green]\n"
                )

                print("Write a vector like [1, 2, 3], then perform operations!")
                print(
                    "\n- Use [bright_magenta]@[/bright_magenta] to calculate the dot product"
                )
                print("- cross(smth, smth) to calculate the cross product")
                print(
                    "- A pair of [bright_magenta]< >[/bright_magenta]s encasing a vector to calculate it's norm!"
                )
                print('Enter any expression to start ("q" to quit):\n')

                while True:
                    expr = input(key_bindings=key_bindings)
                    try:
                        if expr != "q":
                            expr = (
                                (
                                    (
                                        (expr.replace("[", "array([")).replace(
                                            "]", "])"
                                        )
                                    ).replace("<", "norm(")
                                ).replace(">", ")")
                            ).strip(" ")
                            result = eval(replace_graph(expr))
                            print("\n[bright_yellow]Result: [/bright_yellow]", end="")
                            pprint(result)
                            print()

                        else:
                            break
                    except:
                        logging.error(f'Could not parse: "{expr}"\n')

                print()

            elif cmd == "3":
                nprint(
                    "\n[bold bright_green](Current screen: Matrix Calculation Screen)[/bold bright_green]\n"
                )

                print("Write a matrix like [[1, 2], [3, 4]], then perform operations!")
                print(
                    "- Use [bright_magenta]@[/bright_magenta] to calculate the dot product"
                )
                print(
                    "- A pair of [bright_magenta]| |[/bright_magenta] to calculate the determinant."
                )
                print('Enter any expression to start ("q" to quit):\n')

                while True:
                    expr = input(key_bindings=key_bindings)
                    try:
                        if expr != "q":
                            expr = (
                                (
                                    (
                                        (expr.replace("[[", "array([[")).replace(
                                            "]]", "]])"
                                        )
                                    ).replace("|a", "det(a")
                                ).replace(")|", "))")
                            ).strip(" ")
                            result = eval(replace_graph(expr))
                            print("\n[bright_yellow]Result: [/bright_yellow]\n")
                            pprint(result)
                            print()

                        else:
                            break
                    except:
                        print("\n")
                        logging.error(f'Could not parse: "{expr}"\n')

            elif cmd == "4":
                print("\n[bright_yellow]Exiting AlgCalc ... ... ...[/bright_yellow]")
                break

            else:
                print("\n")
                logging.warning(f'Invalid command: "{cmd}"')

        except:
            print("\n")
            logging.error("UnknownError: An unknown error occured.")
            nprint("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def settings():
    while True:
        settings_path = (
            os.path.dirname(os.path.abspath(__file__)) + "/texts/settings.txt"
        )
        settings = open(settings_path, mode="r")

        for line in settings.readlines():
            line = line.rstrip()
            if ("---" in line) or ("|" in line):
                nprint("[gold1]" + line + " [/gold1]")
            else:
                nprint(line)
        print("\n[bold green](Current Screen: Settings Screen)[/bold green]\n")
        print("[bright_purple]Enter Command: [/bright_purple]", end="")
        cmd = input()

        if cmd == "print":
            print(
                "\n[bold bright_green](Current Screen: Print Settings Screen)[/bold bright_green]\n"
            )

            global print_option
            print_option = input(
                'Set print mode: "p" (Sympy Pretty Print) or "n" (Normal Print): '
            )
            print(f'\nPrinting mode set to: "{print_option}"')

        elif cmd == "graph":
            print(
                "\n[bold bright_green](Current Screen: Graph Settings Screen)[/bold bright_green]\n"
            )

            global graph_option
            graph_option = input(
                'Set graph mode: "f" (Fixed graph view) or "a" (Adjusted graph view): '
            )
            print(f'\nGraph mode set to: "{graph_option}"')

        elif cmd == "date":
            print(
                "\n[bold bright_green](Current Screen: Date Settings Screen)[/bold bright_green]\n"
            )

            global date_option
            date_option = input(
                'Set date mode: "1" (YY/MM/DD) or "2" (YY/MM/DD/HH/MM/SS): '
            )
            print(f'\nDate mode set to: "{date_option}"')

        elif cmd == "style":
            print(
                "\n[bold bright_green](Current Screen: Graph Style Settings Screen)[/bold bright_green]\n"
            )

            print("You currently have these available styles:\n")

            style_list = plt.style.available
            styles = ", ".join(style_list)
            nprint("[bright_yellow]" + styles + " [/bright_yellow]")

            style = input('\nChoose a style to apply (type "default" to reset): ')

            plt.style.use(style)

            print(f'Graph style set to "{style}".')

        elif cmd == "exit":
            print("\n[bright_yellow]Exiting settings ... ... ...[/bright_yellow]")
            break

        else:
            print("\n")
            logging.warning(f'Invalid command:"{cmd}"')


main()
# You've reached the end of the file!
