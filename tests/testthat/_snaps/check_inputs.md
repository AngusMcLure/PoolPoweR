# check_geq()

    Code
      check_geq("pool_size", "chr")
    Message <cliMessage>
      i pool_size needs to be a numeric value 0 or greater.
      x "chr" is a character.

---

    Code
      check_geq("pool_size", -1)
    Message <cliMessage>
      i pool_size needs to be a numeric value 0 or greater.
      x -1 is < 0

---

    Code
      check_geq("max_s", 0)
    Message <cliMessage>
      i max_s needs to be a numeric value 1 or greater.
      x 0 is < 1

