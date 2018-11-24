from math import floor, ceil
from sys import stdout as so
from bisect import bisect


def encode(x, p):

    precision = 32
    one = int(2**precision - 1)
    quarter = int(ceil(one/4))
    half = 2*quarter
    threequarters = 3*quarter

    p = dict([(a, p[a]) for a in p if p[a] > 0])

    # Compute cumulative probability as in Shannon-Fano but WITHOUT
    # sorting the probability distrbution
    # ...
    # ...
    # ...
    # ...
    f = dict([(a, mf) for a, mf in zip(p, f)])

    y = []  # initialise output list
    lo, hi = 0, one  # initialise lo and hi to be [0,1.0)
    straddle = 0  # initialise the straddle counter to 0

    for k in range(len(x)):  # for every symbol

        # arithmetic coding is slower than vl_encode, so we display a "progress bar"
        # to let the user know that we are processing the file and haven't crashed...
        if k % 100 == 0:
            so.write('Arithmetic encoded %d%%    \r' % int(floor(k/len(x)*100)))
            so.flush()

        # 1) calculate the interval range to be the difference between hi and lo and
        # add 1 to the difference. The added 1 is necessary to avoid rounding issues
        # lohi_range = ....

        # 2) narrow the interval end-points [lo,hi) to the new range [f,f+p]
        # within the old interval [lo,hi], being careful to round 'innwards' so
        # the code remains prefix-free (you want to use the functions ceil and
        # floor). This will require two instructions. Note that we start computing
        # the new 'lo', then compute the new 'hi' using the scaled probability as
        # the offset from the new 'lo' to the new 'hi'
        # ...
        # ...

        if (lo == hi):
            raise NameError('Zero interval!')

        # Now we need to re-scale the interval if its end-points have bits in common,
        # and output the corresponding bits where appropriate. We will do this with an
        # infinite loop, that will break when none of the conditions for output / straddle
        # are fulfilled
        while True:
            if hi < half:  # if lo < hi < 1/2
                # stretch the interval by 2 and output a 0 followed by 'straddle' ones (if any)
                # and zero the straddle after that. In fact, HOLD OFF on doing the stretching:
                # we will do the stretching at the end of the if statement
                # ...  append a zero to the output list y
                # ...  extend by a sequence of 'straddle' ones
                # ...  zero the straddle counter
            elif lo >= half:  # if hi > lo >= 1/2
                # stretch the interval by 2 and substract 1, and output a 1 followed by 'straddle'
                # zeros (if any) and zero straddle after that. Again, HOLD OFF on doing the stretching
                # as this will be done after the if statement, but note that 2*interval - 1 is equivalent
                # to 2*(interval - 1/2), so for now just substract 1/2 from the interval upper and lower
                # bound (and don't forget that when we say "1/2" we mean the integer "half" we defined
                # above: this is an integer arithmetic implementation!
                # ...  append a 1 to the output list y
                # ...  extend 'straddle' zeros
                # ...  reset the straddle counter
                # ...
                # ...  substract half from lo and hi
            elif lo >= quarter and hi < threequarters:  # if 1/4 < lo < hi < 3/4
                # we can increment the straddle counter and stretch the interval around
                # the half way point. This can be impemented again as 2*(interval - 1/4),
                # and as we will stretch by 2 after the if statement all that needs doing
                # for now is to subtract 1/4 from the upper and lower bound
                # ...  increment straddle
                # ...
                # ...  subtract 'quarter' from lo and hi
            else:
                break  # we break the infinite loop if the interval has reached an un-stretchable state
            # now we can stretch the interval (for all 3 conditions above) by multiplying by 2
            # ...  multiply lo by 2
            # ...  multiply hi by 2 and add 1 (I DON'T KNOW WHY +1 IS NECESSARY BUT IT IS. THIS IS MAGIC.
            #      A BOX OF CHOCOLATES FOR ANYONE WHO GIVES ME A WELL ARGUED REASON FOR THIS... It seems
            #      to solve a minor precision problem.)

    # termination bits
    # after processing all input symbols, flush any bits still in the 'straddle' pipeline
    straddle += 1  # adding 1 to straddle for "good measure" (ensures prefix-freeness)
    if lo < quarter:  # the position of lo determines the dyadic interval that fits
        # ...  output a zero followed by "straddle" ones
        # ...
    else:
        # ...  output a 1 followed by "straddle" zeros

    return(y)


def decode(y, p, n):
    precision = 32
    one = int(2**precision - 1)
    quarter = int(ceil(one/4))
    half = 2*quarter
    threequarters = 3*quarter

    p = dict([(a, p[a]) for a in p if p[a] > 0])

    alphabet = list(p)
    f = [0]
    for a in p:
        f.append(f[-1]+p[a])
    f.pop()

    p = list(p.values())

    y.extend(precision*[0])  # dummy zeros to prevent index out of bound errors
    x = n*[0]  # initialise all zeros

    # initialise by taking first 'precision' bits from y and converting to a number
    value = int(''.join(str(a) for a in y[0:precision]), 2)
    y_position = precision  # position where currently reading y
    lo, hi = 0, one

    x_position = 0
    while 1:
        if x_position % 100 == 0:
            so.write('Arithmetic decoded %d%%    \r' % int(floor(x_position/n*100)))
            so.flush()

        lohi_range = hi - lo + 1
        a = bisect(f, (value-lo)/lohi_range) - 1
        x[x_position] = alphabet[a]

        lo = lo + int(ceil(f[a]*lohi_range))
        hi = lo + int(floor(p[a]*lohi_range))
        if (lo == hi):
            raise NameError('Zero interval!')

        while True:
            if hi < half:
                # do nothing
                pass
            elif lo >= half:
                lo = lo - half
                hi = hi - half
                value = value - half
            elif lo >= quarter and hi < threequarters:
                lo = lo - quarter
                hi = hi - quarter
                value = value - quarter
            else:
                break
            lo = 2*lo
            hi = 2*hi + 1
            value = 2*value + y[y_position]
            y_position += 1
            if y_position == len(y):
                break

        x_position += 1
        if x_position == n or y_position == len(y):
            break

    return(x)
