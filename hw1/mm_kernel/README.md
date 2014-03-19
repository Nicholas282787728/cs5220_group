# Matrix multiply kernel tuning

There are fundamentally two ways to lose performance for matrix
multiplication, mirroring the two points mentioned in lecture 3:

1. We can lose performance because of poor cache utilization.
2. We can lose performance because of poor processor utilization.

The performance loss due to memory access times may be the more
counterintuitive of these two issues, but it is also the one that is
perhaps most straightforward to deal with by techniques such as
blocking and copy optimization for improved layout and alignment.  The
performance loss due to poor processor utilization, though, is a bit
more counterintuitive.

## Building on fast kernels

It frequently does not make sense to expend the effort required to get
the best possible performance.  Pretty good performance may be good enough,
and there's often a law of diminishing returns -- once you've eliminated the
easy bottlenecks, the only opportunities that remain require a great deal
of work.  If we're going to make our code a mess in order to squeeze out as
much performance as possible, it makes good sense to at least try to contain
the mess.  That suggests that if it's necessary to write black magic code,
we should at the very least try to tuck it into a small kernel so that
a later reader does not have to exhaust himself.

If we're going to start from a fast kernel, it also makes sense to
give ourselves the instruments to see whether we're succeeding in
creating something fast from which we can build higher-level codes.
That means we need a timer.  The `ktimer.c` module is exactly such a
timer.  We assume a fixed-size kernel matrix multiply (`kdgemm`) that
takes in pre-digested data that is aligned and laid out however we'd
like.  The kernel routine has to tell us the size of the matrices it
will accept, and from there we run a timing loop to compute an
effective gigaflop rate.

It's also useful to automatically test the kernel in this harness.

## Principles for a fast kernel

There are a few principles involved in getting good instruction-level
parallelism, and then there is some black magic.  Some principles are:

1. Ask for the compiler's help!  At the very least, one should use an
   appropriate set of base optimizations (`-O2` or `-O3`) and tell the
   compiler what type of instruction set it should target for tuning.
   Note that `-O3` is *not* necessarily faster than `-O2`!

2. Organize for independent subcomputations.  For example, if we
   wanted to take the sum of four numbers, the expression
   `(a+b)+(c+d)` involves two adds that can be scheduled in parallel
   by the processor, while in the expression `((a+b)+c)+d`, though
   mathematically equivalent, the processor is forced to serialize the
   additions.

3. Provide the compiler with a good mix of instructions.  On the
   Nehalem, the best we can do is to retire one vector add and one
   vector multiply at each cycle.  That means the compiler needs to
   have a mix of adds and multiplies visible when it is doing
   instruction scheduling.

In an ideal world, a vectorizing compiler would be smart enough to
carry things from here.  We don't live in an ideal world, though, an
while the Intel compiler optimizer does some remarkable things (and
the GNU compiler optimizers have gotten better, too), they still can
miss some things.  Using the SSE or AVX primitives directly can boost
performance significantly over the default vectorization.
