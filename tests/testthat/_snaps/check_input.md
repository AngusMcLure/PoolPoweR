# check_geq()

    Code
      check_geq("pool_size", "chr")
    Message <cliMessage>
      i pool_size must be a numeric value 0 or greater.
      x "chr" is a character.

---

    Code
      check_geq("pool_size", -1)
    Message <cliMessage>
      i pool_size must be a numeric value 0 or greater.
      x -1 is < 0

---

    Code
      check_geq("max_s", 0)
    Message <cliMessage>
      i max_s must be a numeric value 1 or greater.
      x 0 is < 1

# check_in_range()

    Code
      check_in_range("prevalence", -1)
    Message <cliMessage>
      i prevalence must be a numeric value between 0 and 1, inclusive.
      x -1 is < 0

---

    Code
      check_in_range("prevalence", 1.1)
    Message <cliMessage>
      i prevalence must be a numeric value between 0 and 1, inclusive.
      x 1.1 is > 1

# check_rho()

    Code
      check_rho(-1)
    Message <cliMessage>
      i rho must be a numeric value between 0 and 1, or NA
      x -1 is < 0

---

    Code
      check_rho(2)
    Message <cliMessage>
      i rho must be a numeric value between 0 and 1, or NA
      x 2 is > 1

---

    Code
      check_rho("chr")
    Message <cliMessage>
      i rho must be a numeric value between 0 and 1, or NA
      x "chr" is a character.

