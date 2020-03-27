--------------------------------------------------
-- tables
--------------------------------------------------

--
-- Name: organism_id_seq; Type: SEQUENCE; Schema: public; Owner: db_owner
--

CREATE SEQUENCE organism_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.organism_id_seq OWNER TO db_owner;

--
-- Name: organism; Type: TABLE; Schema: public; Owner: db_owner; Tablespace: 
--

CREATE TABLE organism (
    id integer DEFAULT nextval('organism_id_seq'::regclass) NOT NULL,
    taxid integer NOT NULL,
    name character varying(255) NOT NULL,
    dbkey character varying(31),
    genome_size integer,
    genome_type genome_type
);


ALTER TABLE public.organism OWNER TO db_owner;

--
-- Name: prep; Type: TABLE; Schema: public; Owner: db_owner; Tablespace: 
--

CREATE TABLE prep (
    id integer NOT NULL,
    sample_id integer NOT NULL,
    user_id integer NOT NULL,
    prep_date date NOT NULL,
    protocol_id integer NOT NULL,
    primer_pair_id integer NOT NULL,
    index_set character varying(63) NOT NULL,
    library_conc real,
    notes text,
    deleted boolean DEFAULT false
);


ALTER TABLE public.prep OWNER TO db_owner;

--
-- Name: prep_id_seq; Type: SEQUENCE; Schema: public; Owner: db_owner
--

CREATE SEQUENCE prep_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.prep_id_seq OWNER TO db_owner;

--
-- Name: prep_protocol; Type: TABLE; Schema: public; Owner: db_owner; Tablespace: 
--

CREATE TABLE prep_protocol (
    id integer NOT NULL,
    nucleic_acid_isolation nucleic_acid_isolation NOT NULL,
    primary_amplification primary_amplification DEFAULT 'none'::primary_amplification,
    notes character varying(255)
);


ALTER TABLE public.prep_protocol OWNER TO db_owner;

--
-- Name: prep_protocol_id_seq; Type: SEQUENCE; Schema: public; Owner: db_owner
--

CREATE SEQUENCE prep_protocol_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.prep_protocol_id_seq OWNER TO db_owner;

--
-- Name: prep_protocol_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: db_owner
--

ALTER SEQUENCE prep_protocol_id_seq OWNED BY prep_protocol.id;


--
-- Name: prep_run_id_seq; Type: SEQUENCE; Schema: public; Owner: db_owner
--

CREATE SEQUENCE prep_run_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.prep_run_id_seq OWNER TO db_owner;

--
-- Name: prep_run; Type: TABLE; Schema: public; Owner: db_owner; Tablespace: 
--

CREATE TABLE prep_run (
    id integer DEFAULT nextval('prep_run_id_seq'::regclass) NOT NULL,
    prep_id integer NOT NULL,
    flowcell_id character varying(34) NOT NULL,
    raw_reads integer,
    is_paired boolean,
    lane integer
);


ALTER TABLE public.prep_run OWNER TO db_owner;

--
-- Name: sample_id_seq; Type: SEQUENCE; Schema: public; Owner: db_owner
--

CREATE SEQUENCE sample_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.sample_id_seq OWNER TO db_owner;

--
-- Name: sample; Type: TABLE; Schema: public; Owner: db_owner; Tablespace: 
--

CREATE TABLE sample (
    id integer DEFAULT nextval('sample_id_seq'::regclass) NOT NULL,
    parent_id integer,
    name character varying(255) NOT NULL,
    species integer DEFAULT 32644,
    sample_type character varying(50),
    created date,
    user_id integer NOT NULL,
    lab character varying(100),
    project character varying(200),
    experiment_type character varying(200),
    notes character varying(255),
    deleted boolean DEFAULT false,
    host integer DEFAULT 32644,
    budget character(8)
);


ALTER TABLE public.sample OWNER TO db_owner;

--
-- Name: sequencer; Type: TABLE; Schema: public; Owner: db_owner; Tablespace: 
--

CREATE TABLE sequencer (
    id character varying(32) NOT NULL,
    platform character varying(32),
    company character varying(32),
    location character varying(32)
);


ALTER TABLE public.sequencer OWNER TO db_owner;

--
-- Name: sequencing_run; Type: TABLE; Schema: public; Owner: db_owner; Tablespace: 
--

CREATE TABLE sequencing_run (
    id character varying(34) NOT NULL,
    run_date date NOT NULL,
    user_id integer NOT NULL,
    sequencer_id character varying(11) NOT NULL,
    notes text,
    deleted boolean,
    path character varying(63),
    total_reads integer,
    total_preps integer,
    forward_length integer,
    reverse_length integer,
    total_lanes integer
);


ALTER TABLE public.sequencing_run OWNER TO db_owner;

--
-- Name: user_id_seq; Type: SEQUENCE; Schema: public; Owner: db_owner
--

CREATE SEQUENCE user_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.user_id_seq OWNER TO db_owner;

--
-- Name: user; Type: TABLE; Schema: public; Owner: db_owner; Tablespace: 
--

CREATE TABLE user (
    id integer DEFAULT nextval('user_id_seq'::regclass) NOT NULL,
    username character varying(8) NOT NULL,
    shortname character(32),
    name character(128) NOT NULL,
    email character(128),
    active boolean,
    group character(3),
    department character(12)
);


ALTER TABLE public.user OWNER TO db_owner;


--------------------------------------------------
-- foreign keys
--------------------------------------------------


--
-- Name: prep_protocol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY prep
    ADD CONSTRAINT prep_protocol_id_fkey FOREIGN KEY (protocol_id) REFERENCES prep_protocol(id);


--
-- Name: prep_sample_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY prep
    ADD CONSTRAINT prep_sample_id_fkey FOREIGN KEY (sample_id) REFERENCES sample(id);


--
-- Name: prep_run_flowcell_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY prep_run
    ADD CONSTRAINT prep_run_flowcell_id_fkey FOREIGN KEY (flowcell_id) REFERENCES sequencing_run(id);


--
-- Name: prep_run_prep_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY prep_run
    ADD CONSTRAINT prep_run_prep_id_fkey FOREIGN KEY (prep_id) REFERENCES prep(id);


--
-- Name: prep_sample_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY prep_old
    ADD CONSTRAINT prep_sample_id_fkey FOREIGN KEY (sample_id) REFERENCES sample(id);


--
-- Name: prep_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY prep
    ADD CONSTRAINT prep_user_id_fkey FOREIGN KEY (user_id) REFERENCES user(id);


--
-- Name: sample_host_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY sample
    ADD CONSTRAINT sample_host_fkey FOREIGN KEY (host) REFERENCES organism(taxid);


--
-- Name: sample_parent_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY sample
    ADD CONSTRAINT sample_parent_id_fkey FOREIGN KEY (parent_id) REFERENCES sample(id);


--
-- Name: sample_species_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY sample
    ADD CONSTRAINT sample_species_fkey FOREIGN KEY (species) REFERENCES organism(taxid);


--
-- Name: sample_user_id_fkey1; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY sample
    ADD CONSTRAINT sample_user_id_fkey1 FOREIGN KEY (user_id) REFERENCES user(id);


--
-- Name: sequencing_run_sequencer_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY sequencing_run
    ADD CONSTRAINT sequencing_run_sequencer_id_fkey FOREIGN KEY (sequencer_id) REFERENCES sequencer(id);


--
-- Name: sequencing_run_user_id_fkey1; Type: FK CONSTRAINT; Schema: public; Owner: db_owner
--

ALTER TABLE ONLY sequencing_run
    ADD CONSTRAINT sequencing_run_user_id_fkey1 FOREIGN KEY (user_id) REFERENCES user(id);
