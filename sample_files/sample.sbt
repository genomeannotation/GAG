Submit-block ::= {
  contact {
    contact {
      name name {
        last "Geib",
        first "Scott"
      },
      affil std {
        affil "USDA-ARS PBARC",
        div "",
        city "Hilo",
        sub "HI",
        country "USA",
        street "64 Nowelo St",
        email "scott.geib@ars.usda.gov",
        fax "",
        phone "808-959-4335",
        postal-code "96720"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Geib",
            first "Scott",
            initials "S.",
            suffix ""
          }
        }
      },
      affil std {
        affil "USDA-ARS PBARC",
        div "",
        city "Hilo",
        sub "HI",
        country "USA",
        street "64 Nowelo St",
        postal-code "96720"
      }
    }
  },
  subtype new
}

Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Geib",
              first "Scott",
              initials "S.",
              suffix ""
            }
          }
        },
        affil std {
          affil "USDA-ARS PBARC",
          div "",
          city "Hilo",
          sub "HI",
          country "USA",
          street "64 Nowelo St",
          postal-code "96720"
        }
      },
      title "Whole Genome Shotgun Sequencing of Bactrocera dorsalis"
    }
  }
}

Seqdesc ::= user {
  type str "DBLink",
  data {
    {
      label str "BioProject",
      num 1,
      data strs {
        "PRJNA208413"
      }
    },
    {
      label str "BioSample",
      num 1,
      data strs {
        "SAMN02203716"
      }
    }
  }
}

