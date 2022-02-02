{{/* vim: set filetype=mustache: */}}
{{/*
Expand the name of the chart.
*/}}
{{- define "banzai-nres.name" -}}
{{- default .Chart.Name .Values.nameOverride | trunc 63 | trimSuffix "-" -}}
{{- end -}}

{{/*
Create a default fully qualified app name.
We truncate at 63 chars because some Kubernetes name fields are limited to this (by the DNS naming spec).
If release name contains chart name it will be used as a full name.
*/}}
{{- define "banzai-nres.fullname" -}}
{{- if .Values.fullnameOverride -}}
{{- .Values.fullnameOverride | trunc 63 | trimSuffix "-" -}}
{{- else -}}
{{- $name := default .Chart.Name .Values.nameOverride -}}
{{- if contains $name .Release.Name -}}
{{- .Release.Name | trunc 63 | trimSuffix "-" -}}
{{- else -}}
{{- printf "%s-%s" .Release.Name $name | trunc 63 | trimSuffix "-" -}}
{{- end -}}
{{- end -}}
{{- end -}}

{{/*
Create chart name and version as used by the chart label.
*/}}
{{- define "banzai-nres.chart" -}}
{{- printf "%s-%s" .Chart.Name .Chart.Version | replace "+" "_" | trunc 63 | trimSuffix "-" -}}
{{- end -}}

{{/*
Common labels
*/}}
{{- define "banzai-nres.labels" -}}
app.kubernetes.io/name: {{ include "banzai-nres.name" . }}
helm.sh/chart: {{ include "banzai-nres.chart" . }}
app.kubernetes.io/instance: {{ .Release.Name }}
{{- if .Chart.AppVersion }}
app.kubernetes.io/version: {{ .Chart.AppVersion | quote }}
{{- end }}
app.kubernetes.io/managed-by: {{ .Release.Service }}
{{- end -}}


{{/*
Generate the PostgreSQL DB hostname
*/}}
{{- define "banzai-nres.dbhost" -}}
{{- if .Values.postgresql.fullnameOverride -}}
{{- .Values.postgresql.fullnameOverride | trunc 63 | trimSuffix "-" -}}
{{- else if .Values.useDockerizedDatabase -}}
{{- printf "%s-postgresql" .Release.Name -}}
{{- else -}}
{{- required "`postgresql.hostname` must be set when `useDockerizedDatabase` is `false`" .Values.postgresql.hostname -}}
{{- end -}}
{{- end -}}

{{/*
Generate the RabbitMQ queue hostname
*/}}
{{- define "banzai-nres.rabbitmq" -}}
{{- if .Values.rabbitmq.fullnameOverride -}}
{{- .Values.rabbitmq.fullnameOverride | trunc 63 | trimSuffix "-" -}}
{{- else if .Values.useDockerizedRabbitMQ -}}
{{- printf "%s-rabbitmq" .Release.Name -}}
{{- else -}}
{{- required "`rabbitmq.hostname` must be set when `useDockerizedRabbitMQ` is `false`" .Values.rabbitmq.hostname -}}
{{- end -}}
{{- end -}}

{{/*
Define shared environment variables
*/}}
{{- define "banzai-nres.Env" -}}
{{/*
External service URLS
*/ -}}
- name: CONFIGDB_URL
  value: {{ .Values.banzaiNres.configdbUrl | quote }}
- name: OBSERVATION_PORTAL_URL
  value: {{ .Values.banzaiNres.observationPortalUrl | quote }}
{{/*
Ingester environment variables
*/ -}}
- name: API_ROOT
  value: {{ .Values.ingester.apiRoot | quote }}
- name: BUCKET
  value: {{ .Values.ingester.s3Bucket | quote }}
- name: FILESTORE_TYPE
  value: {{ .Values.ingester.filestoreType | quote }}
{{- if .Values.ingester.noMetrics }}
- name: OPENTSDB_PYTHON_METRICS_TEST_MODE
  value: "1"
{{- else }}
- name: INGESTER_PROCESS_NAME
  value: {{ required "ingester.ingesterProcessName must be defined if noMetrics is false" .Values.ingester.ingesterProcessName | quote }}
- name: OPENTSDB_HOSTNAME
  value: {{ required "ingester.opentsdbHostname must be defined if noMetrics is false" .Values.ingester.opentsdbHostname | quote }}
- name: OPENTSDB_PORT
  value: {{ required "ingester.opentsdbPort must be defined if noMetrics is false" .Values.ingester.opentsdbPort | quote }}
{{- end }}
- name: POSTPROCESS_FILES
{{- if .Values.ingester.postProcessFiles }}
  value: "True"
{{- else }}
  value: "False"
{{- end }}
- name: AUTH_TOKEN
  valueFrom:
    secretKeyRef:
      name: banzai-nres-secrets
      key: AUTH_TOKEN
- name: AWS_ACCESS_KEY_ID
  valueFrom:
    secretKeyRef:
      name: banzai-nres-secrets
      key: AWS_ACCESS_KEY_ID
- name: AWS_SECRET_ACCESS_KEY
  valueFrom:
    secretKeyRef:
      name: banzai-nres-secrets
      key: AWS_SECRET_ACCESS_KEY
{{- /*
Optional raw data source environment variables
*/}}
{{- if .Values.banzaiNres.useDifferentArchiveSources }}
- name: RAW_DATA_API_ROOT
  value: {{ required "banzaiNres.rawDataApiRoot is required if banzaiNres.useDifferentArchiveSources is specified." .Values.banzaiNres.rawDataApiRoot | quote }}
- name: RAW_DATA_AUTH_TOKEN
  valueFrom:
    secretKeyRef:
      name: banzai-nres-secrets
      key: RAW_DATA_AUTH_TOKEN
{{- end -}}
{{/*
BANZAI-NRES DB Configuration
*/}}
- name: DB_HOST
  value: {{ include "banzai-nres.dbhost" . | quote }}
- name: DB_PASSWORD
  valueFrom:
    secretKeyRef:
      name: banzai-nres-secrets
      key: postgresqlPassword
- name: DB_USER
  value: {{ .Values.postgresql.postgresqlUsername | quote }}
- name: DB_NAME
  value: {{ .Values.postgresql.postgresqlDatabase | quote }}
- name: DB_ADDRESS
  value: postgresql://$(DB_USER):$(DB_PASSWORD)@$(DB_HOST)/$(DB_NAME)
{{/*
Celery task queue configuration
*/ -}}
- name: RABBITMQ_HOST
  value: {{ include "banzai-nres.rabbitmq" . | quote }}
- name: RABBITMQ_PASSWORD
  valueFrom:
    secretKeyRef:
      name: banzai-nres-secrets
      key: rabbitmq-password
- name: TASK_HOST
  value: amqp://{{ .Values.rabbitmq.rabbitmq.username }}:$(RABBITMQ_PASSWORD)@$(RABBITMQ_HOST)/{{ .Values.rabbitmq.vhost }}
- name: RETRY_DELAY
  value: "600000"
- name: CALIBRATE_PROPOSAL_ID
  value: {{ .Values.banzaiNres.calibrateProposalId | quote }}
- name: FITS_BROKER
  value: {{ .Values.banzaiNres.fitsBroker | quote }}
- name: FITS_EXCHANGE
  value: {{ .Values.banzaiNres.fitsExchange | quote }}
- name: QUEUE_NAME
  value: {{ .Values.banzaiNres.queueName | quote }}
- name: CELERY_TASK_QUEUE_NAME
  value: {{ .Values.banzaiNres.celeryTaskQueueName | quote }}
- name: BANZAI_WORKER_LOGLEVEL
  value: {{ .Values.banzaiNres.banzaiWorkerLogLevel | quote }}
- name: PHOENIX_FILE_LOCATION
  value: {{ .Values.banzaiNres.phoenixFileLocation | quote }}
- name: PHOENIX_MODEL_AWS_ACCESS_KEY_ID
  valueFrom:
    secretKeyRef:
      name: banzai-nres-secrets
      key: PHOENIX_MODEL_AWS_ACCESS_KEY_ID
- name: PHOENIX_MODEL_AWS_SECRET_ACCESS_KEY
  valueFrom:
    secretKeyRef:
      name: banzai-nres-secrets
      key: PHOENIX_MODEL_AWS_SECRET_ACCESS_KEY
{{- end -}}
