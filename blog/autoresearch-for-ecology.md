# Autoresearch for Ecology: What Happens When AI Agents Build Mechanistic Models

*Willem Weertman & Anton Star — March 9, 2026*

---

Andrej Karpathy just released [autoresearch](https://github.com/karpathy/autoresearch), a 630-line Python tool that lets AI agents autonomously iterate on ML experiments. The loop is elegant: human writes research instructions, agent modifies training code, runs a 5-minute experiment, evaluates the result, commits if it improves, repeats. Tobi Lutke ran it overnight and woke up to a 19% improvement in model scores.

We've been running essentially the same loop — but for ecology, not ML. For the past month, an AI research assistant and I have been building a coupled eco-evolutionary epidemiological model for the sunflower sea star (*Pycnopodia helianthoides*) die-off. 896 sites spanning Alaska to Baja California. 5,000 individuals per site. Disease dynamics, host genetics, pathogen evolution, larval dispersal, environmental forcing from real satellite data. The works.

Our iteration loop looks like this: I write biological constraints and structural decisions. The AI proposes parameter configurations, writes calibration configs, launches 24 parallel runs on a 64-core Xeon server (~5 hours each), evaluates RMSE against observed recovery data, and designs the next batch. We've done 154 calibration iterations so far, systematically exploring a parameter space that would have taken a single researcher months to navigate.

**The parallels to autoresearch are obvious. The differences are instructive.**

## Where AI Agents Excel

Most of what got us from iteration W01 to W154 was parameter optimization — sweeping K_half values, tuning environmental pathogen decay rates, finding the right connectivity exponent. This is squarely in autoresearch territory. The AI is tireless, systematic, and keeps meticulous records. It proposes reasonable parameter ranges from literature, runs sensitivity analyses, and maintains a detailed log of every decision. 154 iterations in a month, 24 runs per batch, each tracked to its config file and result. No human could sustain that pace.

The AI also handles an enormous amount of grunt work that would otherwise slow the science to a crawl: writing simulation code, generating publication-quality figures, compiling LaTeX reports, managing a 4,500-line codebase with 922 tests, profiling performance bottlenecks, even reading papers and extracting parameter estimates. These aren't trivial tasks — they're the 80% of research that isn't insight.

## Where the Human Expert Is Irreplaceable (For Now)

But here's what Karpathy's autoresearch doesn't face: *structural ambiguity*. In ML training, the objective is clear (lower loss), the search space is well-defined (architecture, hyperparameters), and the evaluation is immediate (5 minutes). The agent never needs to ask "should I be optimizing this metric at all?"

In mechanistic ecological modeling, the hardest decisions aren't about parameters. They're about *mechanisms*:

- Should we model pathogen evolution as individual strain tracking or community-level adaptation? (We chose community — it's more biologically defensible for environmental bacteria and computationally tractable.)

- Is the Alaska population recovery driven by La Niña cooling events or larval import pulses? (The AI initially built a whole narrative around La Niña. I caught that it was wrong — the 2022 bump is recruitment, not temperature relief. That mistake would have propagated into months of misdirected calibration.)

- Do enclosed waterways trap pathogen (bad for sea stars) or retain larvae (good for recovery)? Both. The balance determines whether fjords are refugia or death traps. This required integrating physical oceanography, larval biology, and Vibrio ecology — knowledge spread across papers that don't cite each other.

- Should we add a salinity mechanism based on the latest fjord refuge paper? Not yet — not before we've published with the current model. Scope creep kills projects.

Each of these decisions saved us weeks of dead-end exploration. They required synthesizing incomplete information across disciplines, applying field intuition, and making judgment calls under uncertainty. That's the hard part of science.

## Three Tiers of Automation

I think there's a useful framework for thinking about where AI agents can operate independently versus where they need human oversight:

**Tier 1: Parameter optimization.** Agents alone, today. This is autoresearch. Given a fixed model structure and a clear objective, search the parameter space systematically. Our calibration sweeps, Karpathy's hyperparameter tuning — same thing, different timescales.

**Tier 2: Mechanism selection in well-characterized systems.** Agents probably can, soon. For a well-studied system — influenza epidemiology, Lotka-Volterra dynamics, standard predator-prey — the literature effectively serves as the expert. An agent that can read 200 papers, extract functional forms, and select from established mechanisms could plausibly build and calibrate a mechanistic model without much human guidance. The search space is mapped. The biology is understood. The papers cite each other.

**Tier 3: Novel mechanism design in poorly-understood systems.** Still needs a human. Our model sits here. We're inventing mechanisms for how *Pycnopodia* pathogen communities evolve thermal tolerance, how wavefront disease spread interacts with VBNC dormancy states, how Allee effects create double extinction thresholds in broadcast spawners. There's no literature review that tells you the answer because the questions haven't been asked before. This requires creative synthesis — the kind of thinking where you connect Vibrio thermal biology to coastal residence times to echinoderm immunology, across papers that have never been in the same bibliography.

## Where This Is Heading

The scaling implication matters: **one domain expert + many AI agents could run many models in parallel.** I don't need to be in every calibration loop. I need to be at the inflection points where the model's structure changes — the moments when a new mechanism is added, a wrong narrative is corrected, or a dead end is recognized. Everything between those decision points can be automated.

Right now, the ratio is roughly one structural decision per week, with dozens of parameter iterations between them. That means ~85% of the research time is already automatable. As AI agents get better at literature synthesis and domain reasoning, Tier 2 will expand, and Tier 3 will shrink. But it won't disappear — because the frontier of science is, by definition, the place where the answers aren't in the papers yet.

Karpathy's autoresearch works because ML training has a fast, clean feedback loop. Ecological modeling has a slow, noisy one — 5 hours per run instead of 5 minutes, and the objective function requires biological judgment to define. But the architecture of the process is the same. And the gains from automation are, if anything, larger in ecology, where the parameter spaces are vast and the human bottleneck is severe.

We're 154 iterations into building a model that may help guide the recovery of a critically endangered species. The AI has written most of the code, run most of the experiments, and generated most of the figures. But the science — the part that decides what the model *should be* — that's still human. For now.

---

*Willem Weertman is a PhD candidate at the University of Washington studying neural systems and behavior, collaborating with the Sea Star Lab at Friday Harbor Laboratories. Anton Star is his AI research assistant. The SSWD-EvoEpi model simulates captive-bred sunflower star reintroduction outcomes for conservation planning.*

*Code: [github.com/anton-openclaw/sswd-evoepi](https://github.com/anton-openclaw/sswd-evoepi)*
