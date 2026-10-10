# aTDS video: narration script for recording

How to record: one audio file per chapter (any format: phone voice memo, .m4a, .mp3, .wav). Read naturally, at your own pace.
Leave about **2 seconds of silence between paragraphs** (between the numbered blocks), so they can be cut apart.
Name the files `ch00.m4a` … `ch10.m4a` (ch00 = the opening). A quiet room, phone about 20 cm from your mouth.
You may change wording if it feels unnatural in your mouth; keep the numbers.

## ch00 · Opening

**1.** Here's something you've felt without noticing. Your heart and your breathing aren't independent. Breathe in, and your heart speeds up a little. But not instantly. It happens a moment later, with a delay.

**2.** Physiologists want to know which systems in the body are linked like this, and when. One popular tool for that is called Time Delay Stability. TDS, for short.

**3.** In this video, we'll look at a version that picks its own settings. We'll call it adaptive TDS. But to get there, we first need to see how TDS itself thinks.

## ch01 · Finding the delay

**4.** Okay. Take a one minute piece of two signals. Is the second one just a delayed copy of the first? Well, let's find out. We slide it in time.

**5.** At every shift, we multiply the two curves point by point, and add it all up. When the shapes line up, most of those products are positive, so the sum gets big. This running score is called the cross correlation.

**6.** And its peak tells us the delay. We call it tau zero. Here, the second signal follows the first by 6 seconds.

## ch02 · Stable delays

**7.** Now, let's do that again. And again. Window after window, along the whole recording. Every window gives us one delay. One dot.

**8.** If two systems really are coupled, the dots line up at the same delay. If they're not, they just jump around.

**9.** TDS turns this into a simple rule. A window counts as stable when at least four out of five neighbouring delays agree, within a small tolerance. In the published method, that's plus or minus one second.

**10.** Let's fill in the stable windows. The TDS score is just the share of windows that are stable.

## ch03 · The hidden assumption

**11.** Now, those settings, a sixty second window, a thirty second step, plus or minus one second, were chosen for sleep recordings sampled once per second. And they quietly assume something. They assume that one minute holds plenty of independent information.

**12.** Look at a fast signal. Every few seconds, it's doing something new. Now look at a slow one. Neighbouring seconds have almost the same value. So, sixty points, but really only a couple of facts.

**13.** And with that little information, the delay in each window is close to random. Here are two slow signals that really are coupled, with a ten second delay. Watch what classic TDS does. The dots scatter. Its score is under three percent, barely above chance.

## ch04 · Step 1: measure memory

**14.** So here's the idea behind adaptive TDS. Before anything else, measure memory. Put a signal next to a copy of itself, shifted by k seconds, and ask: how similar are they?

**15.** With no shift, they match perfectly. A correlation of one. As the shift grows, the match fades. That curve is the autocorrelation. And notice, a fast signal forgets within seconds.

**16.** The Bartlett factor adds up the squared autocorrelation, until the curve first crosses zero. B equals one, plus twice the sum of r squared. Think of it as: how many samples in a row are worth one independent sample?

**17.** For the fast signal, B is about 2. For the slow one, it's about 31. So 31 seconds of that slow signal carry roughly one second's worth of new information.

## ch05 · Step 2: size the window

**18.** Step two. Make the window long enough to hold about thirty independent samples. So the window is thirty times B. And we take B from the slowest signal in the recording, so that every pair is measured with the same ruler.

**19.** For our slow pair, that's a window of about 15 minutes. Not one.

**20.** Of course, a longer window means coarser timing. So you set a cap: the longest window you're willing to accept. For the sleep data, we used five minutes.

## ch06 · Step 3: the delay search

**21.** Step three. The window moves forward by half its length, and the delay search grows with it, up to half the window. Slow systems can have long delays, and a fixed thirty second search would simply never find them.

**22.** There's one exception, and that's rhythms. If a signal repeats every T seconds, then a delay of just over half a period looks exactly like a negative delay of just under half a period.

**23.** See? The two overlays look the same. So the search stops at half of the shortest rhythm in the data.

## ch07 · Step 4: the tolerance

**24.** Step four is the subtle one. The tolerance. With a long window, plus or minus one second is way too strict, because a real delay wanders a little. But if you make the tolerance too loose, chance alone starts to look like stability.

**25.** So, adaptive TDS asks the data. It takes one signal and rotates it in time by a large random amount. Whatever falls off the end wraps around to the start. Each signal keeps its own shape and rhythm. But any real coupling between them is gone.

**26.** Then it runs TDS on these fake pairs, while slowly widening the tolerance. And it keeps the widest tolerance at which the fake pairs still look coupled no more than five percent of the time.

**27.** For our slow pair, that's plus or minus 43 seconds. That's wide. But it's earned. With it, chance stays at or below five percent.

## ch08 · Putting it together

**28.** Now let's put it all together, on that same slow pair. Classic TDS, with one minute windows. The dots scatter. Adaptive TDS, with 15 minute windows. And now, the dots lock onto the true ten second delay. Half of its windows are stable, against a chance level of zero.

**29.** That was just one example. But in earlier tests, with pairs whose coupling and delay were known, the pattern held. As the signals get slower, classic TDS falls to chance. Adaptive TDS stays close to the best possible setting. That's the one you'd pick if you already knew the answer.

## ch09 · Real data

**30.** What about real data? On sleep recordings from Bashan's study, thirty five people, both methods find the known result. The body's network is weakest in deep sleep. Here, adaptive TDS doesn't add much. And that makes sense, because the classic settings were designed for exactly this data.

**31.** And on one of your runs, the fixed app setting found no stable delay between speed and heart rate. Adaptive TDS found 5 stable windows, with heart rate following speed by about 16 seconds.

## ch10 · The price, and a summary

**32.** Of course, nothing is free. A long window blurs short events. A five minute window can't isolate a two minute stretch of being awake. And a wide tolerance blurs the exact delay.

**33.** So, let's recap. Adaptive TDS, in four steps. One: measure each signal's memory. Two: size the window to hold enough independent samples. Three: let the delay search grow with it. Four: calibrate the tolerance on pairs where the coupling has been destroyed. And then, choose the cap that matches the time resolution you need.

**34.** And that's adaptive Time Delay Stability. Thanks for watching.
